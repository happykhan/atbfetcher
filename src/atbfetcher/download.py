"""Download and extract genome assemblies from ATB tarballs or AWS S3.

ATB stores assemblies in tar.xz archives on OSF, and also provides individual
gzipped FASTA files on a public AWS S3 bucket. This module supports both
download methods and can auto-select the faster one based on genome count.
"""

import hashlib
import logging
import os
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

DEFAULT_THREADS = max(1, os.cpu_count() - 1) if os.cpu_count() else 1

# AWS S3 public bucket base URL for individual assemblies
AWS_BASE_URL = "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com"

# Empirical timing constants (seconds), measured on typical connections
_AWS_PER_FILE_SECS = 0.4       # avg time to download one .fa.gz from AWS (~800KB)
_TARBALL_DOWNLOAD_SECS = 60.0  # avg time to download one ~3GB tarball (uncached)
_TARBALL_EXTRACT_SECS = 25.0   # avg time to decompress + scan one tarball


def estimate_download_time(
    n_genomes: int,
    n_tarballs: int,
) -> tuple[str, float, float]:
    """Estimate download time for AWS vs tarball and pick the faster method.

    Assumes uncached tarballs — if the user has cached tarballs they can
    force ``--source osf`` to take advantage of that.

    Parameters
    ----------
    n_genomes : int
        Number of genomes to download.
    n_tarballs : int
        Number of unique tarballs needed.

    Returns
    -------
    tuple[str, float, float]
        (recommended_method, aws_estimate_secs, tarball_estimate_secs)
    """
    aws_time = n_genomes * _AWS_PER_FILE_SECS
    tarball_time = n_tarballs * (_TARBALL_DOWNLOAD_SECS + _TARBALL_EXTRACT_SECS)
    method = "aws" if aws_time <= tarball_time else "osf"
    return method, aws_time, tarball_time


def resolve_tarballs(
    sample_ids: list[str], file_list_df: pd.DataFrame
) -> dict[str, list[dict]]:
    """Map sample IDs to their containing tarballs.

    Groups samples by tarball so each archive is downloaded only once.

    Parameters
    ----------
    sample_ids : list[str]
        Sample accession IDs to look up.
    file_list_df : pd.DataFrame
        File list DataFrame with columns: sample, filename_in_tar_xz,
        tar_xz, tar_xz_url, tar_xz_md5.

    Returns
    -------
    dict
        Mapping of tarball name to list of dicts with sample info:
        ``{tarball: [{sample, filename, url, md5}, ...]}``.
    """
    matched = file_list_df[file_list_df["sample"].isin(sample_ids)]

    if len(matched) < len(sample_ids):
        missing = set(sample_ids) - set(matched["sample"])
        logger.warning(
            "%d samples not found in file list: %s",
            len(missing),
            ", ".join(list(missing)[:5]),
        )

    tarballs: dict[str, list[dict]] = {}
    for _, row in matched.iterrows():
        tarball = row["tar_xz"]
        entry = {
            "sample": row["sample"],
            "filename": row["filename_in_tar_xz"],
            "url": row["tar_xz_url"],
            "md5": row["tar_xz_md5"],
        }
        tarballs.setdefault(tarball, []).append(entry)

    return tarballs


def download_tarball(url: str, dest: Path, md5: str | None = None) -> Path:
    """Download a tarball with a progress bar and optional MD5 verification.

    Parameters
    ----------
    url : str
        URL to download from.
    dest : Path
        Local path to save the file.
    md5 : str, optional
        Expected MD5 hash for verification.

    Returns
    -------
    Path
        Path to the downloaded file.

    Raises
    ------
    ValueError
        If MD5 verification fails.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading %s", url)
    response = requests.get(url, stream=True, timeout=300)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    hasher = hashlib.md5()

    with (
        open(dest, "wb") as f,
        tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            desc=dest.name,
        ) as pbar,
    ):
        for chunk in response.iter_content(chunk_size=65536):
            f.write(chunk)
            hasher.update(chunk)
            pbar.update(len(chunk))

    if md5 and hasher.hexdigest() != md5:
        dest.unlink()
        raise ValueError(
            f"MD5 mismatch for {dest.name}: "
            f"expected {md5}, got {hasher.hexdigest()}"
        )

    return dest


def _decompress_xz(tarball_path: Path, threads: int) -> Path:
    """Decompress .tar.xz to .tar using multi-threaded xz.

    Parameters
    ----------
    tarball_path : Path
        Path to the .tar.xz file.
    threads : int
        Number of threads for xz decompression.

    Returns
    -------
    Path
        Path to the decompressed .tar file.
    """
    tar_path = tarball_path.with_suffix("")  # remove .xz
    if tar_path.exists():
        return tar_path

    logger.info("Decompressing %s with %d threads", tarball_path.name, threads)
    subprocess.run(
        ["xz", "-dk", f"-T{threads}", str(tarball_path)],
        check=True,
        capture_output=True,
    )
    return tar_path


def extract_samples(
    tarball_path: Path,
    sample_filenames: list[str],
    output_dir: Path,
    threads: int = DEFAULT_THREADS,
) -> list[Path]:
    """Extract specific assembly files from a tar.xz archive.

    Uses multi-threaded XZ decompression when threads > 1.

    Parameters
    ----------
    tarball_path : Path
        Path to the tar.xz archive.
    sample_filenames : list[str]
        Filenames within the archive to extract.
    output_dir : Path
        Directory to write extracted files.
    threads : int
        Number of threads for XZ decompression.

    Returns
    -------
    list[Path]
        Paths to extracted assembly files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build lookup: match on both full path and basename
    target_full = set(sample_filenames)
    target_basenames = {Path(f).name for f in sample_filenames}
    extracted = []

    # Use multi-threaded xz decompression if threads > 1 and xz is available
    tar_path = None
    use_subprocess_xz = threads > 1
    if use_subprocess_xz:
        try:
            tar_path = _decompress_xz(tarball_path, threads)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.debug("xz CLI not available, falling back to Python tarfile")
            use_subprocess_xz = False

    open_path = tar_path if use_subprocess_xz else tarball_path
    open_mode = "r:" if use_subprocess_xz else "r:xz"

    with tarfile.open(open_path, open_mode) as tar:
        for member in tar.getmembers():
            basename = Path(member.name).name
            if member.name in target_full or basename in target_basenames:
                member.name = basename
                tar.extract(member, path=output_dir, filter="data")
                extracted.append(output_dir / basename)
                logger.debug("Extracted %s", basename)

    # Clean up decompressed tar
    if tar_path and tar_path.exists():
        tar_path.unlink()

    if len(extracted) < len(sample_filenames):
        found = {p.name for p in extracted}
        expected_basenames = {Path(f).name for f in sample_filenames}
        missing = expected_basenames - found
        logger.warning(
            "%d files not found in tarball: %s",
            len(missing),
            ", ".join(list(missing)[:5]),
        )

    return extracted


def _download_one_aws(sample_id: str, output_dir: Path) -> Path | None:
    """Download a single assembly from AWS S3.

    Parameters
    ----------
    sample_id : str
        Sample accession ID (e.g. SAMN12345678).
    output_dir : Path
        Directory to save the file.

    Returns
    -------
    Path or None
        Path to the downloaded .fa.gz file, or None on failure.
    """
    url = f"{AWS_BASE_URL}/{sample_id}.fa.gz"
    dest = output_dir / f"{sample_id}.fa.gz"

    if dest.exists():
        logger.debug("Already exists: %s", dest.name)
        return dest

    try:
        response = requests.get(url, stream=True, timeout=120)
        response.raise_for_status()

        with open(dest, "wb") as f:
            for chunk in response.iter_content(chunk_size=65536):
                f.write(chunk)

        return dest
    except requests.RequestException as e:
        logger.warning("Failed to download %s from AWS: %s", sample_id, e)
        if dest.exists():
            dest.unlink()
        return None


def fetch_from_aws(
    sample_ids: list[str],
    output_dir: Path,
    max_workers: int = 8,
) -> list[Path]:
    """Download assemblies individually from the ATB AWS S3 bucket.

    Downloads gzipped FASTA files in parallel from the public bucket.

    Parameters
    ----------
    sample_ids : list[str]
        Sample accession IDs to download.
    output_dir : Path
        Directory for downloaded assemblies.
    max_workers : int
        Number of concurrent downloads.

    Returns
    -------
    list[Path]
        Paths to successfully downloaded files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_download_one_aws, sid, output_dir): sid
            for sid in sample_ids
        }
        with tqdm(total=len(sample_ids), desc="AWS downloads", unit="file") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    downloaded.append(result)
                pbar.update(1)

    logger.info("Downloaded %d/%d assemblies from AWS", len(downloaded), len(sample_ids))
    return downloaded


def fetch_assemblies(
    selected_df: pd.DataFrame,
    file_list_df: pd.DataFrame,
    output_dir: Path,
    cache_dir: Path,
    no_cache: bool = False,
    threads: int = DEFAULT_THREADS,
) -> list[Path]:
    """Orchestrate downloading and extracting selected genome assemblies.

    Maps selected samples to tarballs, downloads each tarball (caching
    by default), extracts the requested assemblies, and optionally
    removes tarballs after extraction.

    Parameters
    ----------
    selected_df : pd.DataFrame
        Selected genomes. Must have a ``sample`` column.
    file_list_df : pd.DataFrame
        File list mapping samples to tarballs.
    output_dir : Path
        Directory for extracted assemblies.
    cache_dir : Path
        Directory for caching downloaded tarballs.
    no_cache : bool
        If True, delete tarballs after extraction.
    threads : int
        Number of threads for XZ decompression.

    Returns
    -------
    list[Path]
        Paths to all extracted assembly files.
    """
    output_dir = Path(output_dir)
    cache_dir = Path(cache_dir)
    tarball_cache = cache_dir / "tarballs"
    tarball_cache.mkdir(parents=True, exist_ok=True)

    sample_ids = selected_df["sample"].tolist()
    tarballs = resolve_tarballs(sample_ids, file_list_df)

    all_extracted = []
    total_tarballs = len(tarballs)

    for i, (tarball_name, entries) in enumerate(tarballs.items(), 1):
        logger.info(
            "Processing tarball %d/%d: %s (%d samples)",
            i,
            total_tarballs,
            tarball_name,
            len(entries),
        )

        tarball_path = tarball_cache / tarball_name
        url = entries[0]["url"]
        md5 = entries[0]["md5"]

        # Download if not already cached
        if not tarball_path.exists():
            download_tarball(url, tarball_path, md5)
        else:
            logger.info("Using cached tarball: %s", tarball_path)

        # Extract requested samples
        filenames = [e["filename"] for e in entries]
        extracted = extract_samples(tarball_path, filenames, output_dir, threads)
        all_extracted.extend(extracted)

        # Remove tarball if no_cache
        if no_cache and tarball_path.exists():
            tarball_path.unlink()
            logger.info("Removed tarball: %s", tarball_path)

    logger.info(
        "Fetched %d assemblies from %d tarballs",
        len(all_extracted),
        total_tarballs,
    )

    return all_extracted
