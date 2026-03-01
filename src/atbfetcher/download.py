"""Download and extract genome assemblies from ATB tarballs.

ATB stores assemblies in tar.xz archives on OSF. This module maps sample IDs
to their containing tarballs, downloads them, and extracts only the requested
assembly files. Supports multi-threaded XZ decompression for faster extraction.
"""

import hashlib
import logging
import os
import subprocess
import tarfile
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

DEFAULT_THREADS = max(1, os.cpu_count() - 1) if os.cpu_count() else 1


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
