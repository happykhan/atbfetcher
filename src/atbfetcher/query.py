"""Query the ATB SQLite metadata database.

Allows filtering genomes by species, country, collection year, host,
isolation source, and quality metrics using the comprehensive ATB
metadata SQLite database.
"""

import logging
import sqlite3
import subprocess
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

# OSF download URL for the compressed SQLite database (~2 GB)
SQLITE_URL = "https://osf.io/download/my56u/"
SQLITE_FILENAME = "atb.metadata.202505.sqlite"
SQLITE_COMPRESSED = f"{SQLITE_FILENAME}.xz"


def find_sqlite_db(db_path: Path | None, cache_dir: Path) -> Path | None:
    """Locate the ATB SQLite database.

    Checks the user-provided path first, then the cache directory.

    Parameters
    ----------
    db_path : Path or None
        Explicit path provided by the user.
    cache_dir : Path
        Default cache directory to look in.

    Returns
    -------
    Path or None
        Path to the SQLite database, or None if not found.
    """
    if db_path and db_path.exists():
        return db_path

    default_path = cache_dir / SQLITE_FILENAME
    if default_path.exists():
        return default_path

    return None


def download_sqlite_db(cache_dir: Path) -> Path:
    """Download and decompress the ATB SQLite metadata database.

    Downloads the ~2 GB compressed database from OSF with a progress bar,
    then decompresses it (~27 GB) using the ``xz`` command.

    Parameters
    ----------
    cache_dir : Path
        Directory to store the database.

    Returns
    -------
    Path
        Path to the decompressed SQLite database.

    Raises
    ------
    RuntimeError
        If download fails or ``xz`` decompression fails.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    compressed_path = cache_dir / SQLITE_COMPRESSED
    db_path = cache_dir / SQLITE_FILENAME

    if db_path.exists():
        logger.info("Database already exists at %s", db_path)
        return db_path

    # Download
    logger.info("Downloading ATB metadata database from %s", SQLITE_URL)
    response = requests.get(SQLITE_URL, stream=True, timeout=30)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))

    with (
        open(compressed_path, "wb") as f,
        tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            desc="Downloading database",
        ) as pbar,
    ):
        for chunk in response.iter_content(chunk_size=65536):
            f.write(chunk)
            pbar.update(len(chunk))

    logger.info(
        "Downloaded %s (%.1f GB)", compressed_path.name, compressed_path.stat().st_size / 1e9
    )

    # Decompress
    logger.info("Decompressing %s (this may take a few minutes)...", compressed_path.name)
    try:
        subprocess.run(
            ["xz", "-d", str(compressed_path)],
            check=True,
            capture_output=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "xz is not installed. Install it with your package manager "
            "(e.g. 'brew install xz' or 'apt install xz-utils') and try again, "
            f"or decompress manually: xz -d {compressed_path}"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Failed to decompress database: {exc.stderr.decode()}") from exc

    if not db_path.exists():
        raise RuntimeError(f"Decompression succeeded but {db_path} not found")

    logger.info("Database ready at %s (%.1f GB)", db_path, db_path.stat().st_size / 1e9)
    return db_path


def query_metadata(
    db_path: Path,
    *,
    species: str | None = None,
    country: str | None = None,
    year_from: int | None = None,
    year_to: int | None = None,
    host: str | None = None,
    isolation_source: str | None = None,
    hq_only: bool = True,
    min_completeness: float | None = None,
    max_contamination: float | None = None,
    min_genome_size: int | None = None,
    max_genome_size: int | None = None,
    limit: int | None = None,
    seed: int = 42,
) -> pd.DataFrame:
    """Query the ATB SQLite database with flexible filters.

    Joins assembly, run, ena_202505_used, and optionally checkm2 tables
    to return matching sample accessions with metadata.

    Parameters
    ----------
    db_path : Path
        Path to the ATB SQLite database.
    species : str, optional
        Filter by Sylph species (exact or LIKE match if contains %).
    country : str, optional
        Filter by country (prefix match, e.g. "United Kingdom" matches
        "United Kingdom", "United Kingdom:England", etc.).
    year_from : int, optional
        Minimum collection year (inclusive).
    year_to : int, optional
        Maximum collection year (inclusive).
    host : str, optional
        Filter by host organism (LIKE match).
    isolation_source : str, optional
        Filter by isolation source (LIKE match).
    hq_only : bool
        If True (default), only return genomes passing ATB HQ filter.
    min_completeness : float, optional
        Minimum CheckM2 completeness (Completeness_Specific).
    max_contamination : float, optional
        Maximum CheckM2 contamination.
    min_genome_size : int, optional
        Minimum genome size in bases.
    max_genome_size : int, optional
        Maximum genome size in bases.
    limit : int, optional
        Maximum number of results to return.
    seed : int
        Random seed for reproducible sampling when limit is set.

    Returns
    -------
    pd.DataFrame
        Matching samples with columns: sample, species, hq_filter,
        country, collection_date, host, isolation_source, aws_url,
        and optionally checkm2 metrics.
    """
    conn = sqlite3.connect(str(db_path))

    # Determine whether we need checkm2 join
    need_checkm2 = any(
        [
            min_completeness is not None,
            max_contamination is not None,
            min_genome_size is not None,
            max_genome_size is not None,
        ]
    )

    # Determine whether we need ENA join
    need_ena = any(
        [
            country is not None,
            year_from is not None,
            year_to is not None,
            host is not None,
            isolation_source is not None,
        ]
    )

    # Build SELECT and FROM
    select_cols = [
        "a.sample_accession AS sample",
        "a.sylph_species AS species",
        "a.hq_filter",
        "a.aws_url",
    ]
    from_clause = "FROM assembly a"
    joins = []
    conditions = []
    params = []

    if need_ena:
        joins.append(
            "JOIN run r ON a.sample_accession = r.sample_accession "
            "JOIN ena_202505_used e ON r.run_accession = e.run_accession"
        )
        select_cols.extend(
            [
                "e.country",
                "e.collection_date",
                "e.host",
                "e.isolation_source",
            ]
        )

    if need_checkm2:
        joins.append("JOIN checkm2 c ON a.sample_accession = c.sample_accession")
        select_cols.extend(
            [
                "c.Completeness_Specific",
                "c.Contamination",
                "c.Genome_Size",
                "c.GC_Content",
                "c.Contig_N50",
            ]
        )

    # Assembly must be available on OSF
    conditions.append("a.asm_fasta_on_osf = 1")

    # HQ filter
    if hq_only:
        conditions.append("a.hq_filter = 'PASS'")

    # Species filter
    if species:
        if "%" in species:
            conditions.append("a.sylph_species LIKE ?")
        else:
            # Match exact or with GTDB suffix (_A, _B, etc.)
            conditions.append(
                "(a.sylph_species = ? OR a.sylph_species LIKE ? || '\\_%' ESCAPE '\\')"
            )
            params.append(species)
        params.append(species)

    # Country filter (prefix match for "Country:Region" format)
    if country:
        conditions.append("e.country LIKE ?")
        params.append(f"{country}%")

    # Year range filter (collection_date can be YYYY, YYYY-MM, YYYY-MM-DD)
    if year_from is not None:
        conditions.append("SUBSTR(e.collection_date, 1, 4) >= ?")
        params.append(str(year_from))
    if year_to is not None:
        conditions.append("SUBSTR(e.collection_date, 1, 4) <= ?")
        params.append(str(year_to))

    # Host filter
    if host:
        conditions.append("e.host LIKE ?")
        params.append(f"%{host}%")

    # Isolation source filter
    if isolation_source:
        conditions.append("e.isolation_source LIKE ?")
        params.append(f"%{isolation_source}%")

    # CheckM2 filters
    if min_completeness is not None:
        conditions.append("c.Completeness_Specific >= ?")
        params.append(min_completeness)
    if max_contamination is not None:
        conditions.append("c.Contamination <= ?")
        params.append(max_contamination)
    if min_genome_size is not None:
        conditions.append("c.Genome_Size >= ?")
        params.append(min_genome_size)
    if max_genome_size is not None:
        conditions.append("c.Genome_Size <= ?")
        params.append(max_genome_size)

    # Assemble query
    sql = f"SELECT {', '.join(select_cols)} {from_clause}"
    if joins:
        sql += " " + " ".join(joins)
    if conditions:
        sql += " WHERE " + " AND ".join(conditions)

    # For reproducible random sampling
    if limit is not None:
        sql += " ORDER BY SUBSTR(HEX(RANDOMBLOB(4)), 1, 8)"
        sql += " LIMIT ?"
        params.append(limit)

    logger.debug("SQL: %s", sql)
    logger.debug("Params: %s", params)

    try:
        # Set PRNG seed for reproducible RANDOM() ordering
        if limit is not None:
            conn.execute(f"SELECT RANDOMBLOB({seed})")

        df = pd.read_sql_query(sql, conn, params=params)
    finally:
        conn.close()

    return df


def list_countries(db_path: Path, species: str | None = None) -> list[str]:
    """List available countries in the ENA metadata.

    Parameters
    ----------
    db_path : Path
        Path to the ATB SQLite database.
    species : str, optional
        If provided, only list countries with genomes of this species.

    Returns
    -------
    list[str]
        Sorted list of country names.
    """
    conn = sqlite3.connect(str(db_path))
    try:
        if species:
            sql = (
                "SELECT DISTINCT e.country FROM assembly a "
                "JOIN run r ON a.sample_accession = r.sample_accession "
                "JOIN ena_202505_used e ON r.run_accession = e.run_accession "
                "WHERE a.hq_filter = 'PASS' AND a.sylph_species = ? "
                "AND e.country IS NOT NULL AND e.country != '' "
                "ORDER BY e.country"
            )
            rows = conn.execute(sql, [species]).fetchall()
        else:
            sql = (
                "SELECT DISTINCT e.country FROM ena_202505_used e "
                "WHERE e.country IS NOT NULL AND e.country != '' "
                "ORDER BY e.country"
            )
            rows = conn.execute(sql).fetchall()
    finally:
        conn.close()

    return [r[0] for r in rows]


def list_hosts(db_path: Path, species: str | None = None) -> list[str]:
    """List available host organisms in the ENA metadata.

    Parameters
    ----------
    db_path : Path
        Path to the ATB SQLite database.
    species : str, optional
        If provided, only list hosts with genomes of this species.

    Returns
    -------
    list[str]
        Sorted list of host names.
    """
    conn = sqlite3.connect(str(db_path))
    try:
        if species:
            sql = (
                "SELECT DISTINCT e.host FROM assembly a "
                "JOIN run r ON a.sample_accession = r.sample_accession "
                "JOIN ena_202505_used e ON r.run_accession = e.run_accession "
                "WHERE a.hq_filter = 'PASS' AND a.sylph_species = ? "
                "AND e.host IS NOT NULL AND e.host != '' "
                "ORDER BY e.host"
            )
            rows = conn.execute(sql, [species]).fetchall()
        else:
            sql = (
                "SELECT DISTINCT e.host FROM ena_202505_used e "
                "WHERE e.host IS NOT NULL AND e.host != '' "
                "ORDER BY e.host"
            )
            rows = conn.execute(sql).fetchall()
    finally:
        conn.close()

    return [r[0] for r in rows]
