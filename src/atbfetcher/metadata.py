"""Download, cache, and load ATB metadata files.

Metadata files (species calls, CheckM2 quality, file lists) are downloaded
from OSF on first use and cached locally as Parquet files for fast reloading.
"""

import gzip
import hashlib
import io
import logging
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

# OSF download URLs for ATB metadata files
URLS = {
    "file_list": "https://osf.io/download/69a040c86a4dd653508ac769/",
    "species_calls": "https://osf.io/download/699f15825f818865268accd2/",
    "checkm2": "https://osf.io/download/699f158c37ec474686e2c7d1/",
}

QUALIBACT_URL = "https://static.qualibact.org/static/summary/filtered_metrics.csv"

DEFAULT_CACHE_DIR = Path.home() / ".atbfetcher"


class MetadataCache:
    """Manages downloading and caching of ATB metadata files.

    Files are downloaded as gzipped TSVs from OSF, then stored locally as
    Parquet files for fast subsequent loads.

    Parameters
    ----------
    cache_dir : Path
        Directory to store cached files. Created if it doesn't exist.
    no_cache : bool
        If True, skip caching — download fresh each time.
    refresh : bool
        If True, re-download even if cached files exist.
    """

    def __init__(
        self,
        cache_dir: Path = DEFAULT_CACHE_DIR,
        no_cache: bool = False,
        refresh: bool = False,
    ):
        self.cache_dir = Path(cache_dir)
        self.no_cache = no_cache
        self.refresh = refresh
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _parquet_path(self, name: str) -> Path:
        """Return the local Parquet cache path for a given dataset name."""
        return self.cache_dir / f"{name}.parquet"

    def _download_gzipped_tsv(self, url: str, name: str) -> pd.DataFrame:
        """Download a gzipped TSV from a URL and return as a DataFrame.

        Shows a progress bar during download using tqdm.
        """
        logger.info("Downloading %s from %s", name, url)
        response = requests.get(url, stream=True, timeout=120)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))
        chunks = []
        with tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            desc=f"Downloading {name}",
            disable=total_size == 0,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                chunks.append(chunk)
                pbar.update(len(chunk))

        raw_bytes = b"".join(chunks)
        decompressed = gzip.decompress(raw_bytes)
        return pd.read_csv(io.BytesIO(decompressed), sep="\t", low_memory=False)

    def _load_or_download(self, name: str, url: str) -> pd.DataFrame:
        """Load from Parquet cache, or download and cache if needed."""
        parquet_path = self._parquet_path(name)

        # Use cache if available and not refreshing
        if parquet_path.exists() and not self.refresh and not self.no_cache:
            logger.info("Loading cached %s from %s", name, parquet_path)
            return pd.read_parquet(parquet_path)

        # Download fresh
        df = self._download_gzipped_tsv(url, name)

        # Save to cache (unless no_cache mode)
        if not self.no_cache:
            df.to_parquet(parquet_path, index=False)
            logger.info("Cached %s to %s", name, parquet_path)

        return df

    def load_file_list(self) -> pd.DataFrame:
        """Load the file list mapping samples to tarballs.

        Columns after normalization: sample, sylph_species,
        filename_in_tar_xz, tar_xz, tar_xz_url, tar_xz_md5, tar_xz_size_MB.
        """
        return self._load_or_download("file_list", URLS["file_list"])

    def load_species_calls(self, hq_only: bool = True) -> pd.DataFrame:
        """Load Sylph species calls for all ATB samples.

        Raw columns are ``sample_accession`` and ``sylph_species``, renamed
        to ``sample`` and ``species`` for consistency across the codebase.

        Parameters
        ----------
        hq_only : bool
            If True (default), keep only samples that pass ATB's HQ filter
            (completeness >= 90%, contamination <= 5%, assembly length
            100kbp–15Mbp, <= 2000 contigs, N50 >= 2000).
        """
        df = self._load_or_download("species_calls", URLS["species_calls"])

        # Normalize column names for internal consistency
        df = df.rename(columns={
            "sample_accession": "sample",
            "sylph_species": "species",
        })

        # Keep only high-quality Sylph calls
        if hq_only and "HQ" in df.columns:
            df = df[df["HQ"] == "T"].copy()

        return df

    def load_checkm2(self) -> pd.DataFrame:
        """Load CheckM2 quality metrics for all ATB samples.

        Raw column ``sample_accession`` is renamed to ``sample``.
        Columns include: sample, Completeness_Specific, Contamination,
        Genome_Size, GC_Content, Contig_N50, and more.
        """
        df = self._load_or_download("checkm2", URLS["checkm2"])

        # Normalize sample ID column
        df = df.rename(columns={"sample_accession": "sample"})

        return df


def load_qualibact_cutoffs(url: str = QUALIBACT_URL) -> dict:
    """Download and parse Qualibact per-species quality cutoffs.

    Returns a nested dict: {species: {metric: (lower_bound, upper_bound)}}.

    Parameters
    ----------
    url : str
        URL to the Qualibact filtered_metrics.csv file.

    Returns
    -------
    dict
        Mapping of species name to dict of metric bounds.
    """
    logger.info("Downloading Qualibact cutoffs from %s", url)
    response = requests.get(url, timeout=60)
    response.raise_for_status()

    df = pd.read_csv(io.StringIO(response.text))

    cutoffs: dict = {}
    for _, row in df.iterrows():
        species = row["species"]
        metric = row["metric"]
        lower = row["lower_bounds"]
        upper = row["upper_bounds"]

        if species not in cutoffs:
            cutoffs[species] = {}
        cutoffs[species][metric] = (lower, upper)

    return cutoffs
