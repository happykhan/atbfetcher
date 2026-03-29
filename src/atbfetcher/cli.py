"""Command-line interface for atbfetcher.

Provides subcommands to fetch genomes from AllTheBacteria:
- ``species``: Subsample by species with stratified sampling
- ``mlst``: Select by MLST sequence types
- ``accessions``: Fetch specific accessions
- ``list-species``: Print available species names
- ``species-count``: Show genome counts per species
- ``query``: Flexible metadata query with optional download
- ``info``: Per-sample detail view
- ``summarise``: Group-by counts and stats
- ``mlst-query``: Filter samples by ST / scheme / species (no download)
- ``config``: Manage user configuration
"""

import difflib
import logging
import sqlite3
import sys
import tomllib
from pathlib import Path

import click
import pandas as pd
from rich.logging import RichHandler

from atbfetcher.config import CONFIG_PATH, default_config, load_config, write_config
from atbfetcher.download import (
    DEFAULT_THREADS,
    estimate_download_time,
    fetch_assemblies,
    fetch_from_aws,
    resolve_tarballs,
)
from atbfetcher.metadata import DEFAULT_CACHE_DIR, MetadataCache, load_qualibact_cutoffs
from atbfetcher.mlst import STRATEGIES, filter_by_mlst, load_suspect_contaminations
from atbfetcher.output import FORMATS, print_dataframe, resolve_format
from atbfetcher.plotting import plot_selection
from atbfetcher.quality import filter_by_quality
from atbfetcher.query import (
    SQLITE_FILENAME,
    SQLITE_URL,
    download_sqlite_db,
    find_sqlite_db,
    list_countries,
    list_hosts,
    query_metadata,
)
from atbfetcher.sampling import stratified_sample
from atbfetcher.species import (
    clean_species_name,
    get_samples_for_species,
    is_placeholder_species,
    list_species,
)


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging with rich for colourful output."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%H:%M:%S]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
    )


# -- Summary table --

# Assembly stat columns to include in the summary TSV
_ASSEMBLY_STATS = [
    "sample",
    "species",
    "Completeness_Specific",
    "Contamination",
    "Genome_Size",
    "GC_Content",
    "Contig_N50",
]

# Extra MLST columns when available
_MLST_COLS = ["mlst_scheme", "mlst_st", "mlst_status"]


def _write_summary(selected_df: pd.DataFrame, output_dir: Path, species_name: str) -> Path:
    """Write a TSV summary of selected genomes to the output directory."""
    cols = [c for c in _ASSEMBLY_STATS + _MLST_COLS if c in selected_df.columns]
    summary = selected_df[cols].copy()
    safe_name = species_name.replace(" ", "_")
    summary_path = output_dir / f"{safe_name}_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    click.echo(f"  Summary table saved to {summary_path}")
    return summary_path


# -- Shared CLI options --


def cache_options(func):
    """Decorator adding shared cache-related CLI options."""
    func = click.option(
        "--cache-dir",
        type=click.Path(path_type=Path),
        default=DEFAULT_CACHE_DIR,
        show_default=True,
        help="Directory for cached metadata and tarballs.",
    )(func)
    func = click.option(
        "--no-cache",
        is_flag=True,
        default=False,
        help="Skip caching — download fresh data each time.",
    )(func)
    func = click.option(
        "--refresh",
        is_flag=True,
        default=False,
        help="Force re-download of cached metadata.",
    )(func)
    return func


def verbose_option(func):
    """Decorator adding a --verbose flag."""
    return click.option(
        "--verbose",
        "-v",
        is_flag=True,
        default=False,
        help="Enable verbose (debug) logging.",
    )(func)


def threads_option(func):
    """Decorator adding a --threads option."""
    return click.option(
        "--threads",
        "-t",
        default=DEFAULT_THREADS,
        show_default=True,
        help="Number of threads for XZ decompression.",
    )(func)


def quality_filter_option(func):
    """Decorator adding a --quality-filter option."""
    return click.option(
        "--quality-filter",
        "-q",
        type=click.Choice(["atb", "qualibact", "none"], case_sensitive=False),
        default="atb",
        show_default=True,
        help=(
            "Quality filter: 'atb' uses ATB HQ flag "
            "(completeness>=90%%, contamination<=5%%, etc.), "
            "'qualibact' adds per-species cutoffs from Qualibact on top of ATB HQ, "
            "'none' skips all quality filtering."
        ),
    )(func)


def source_option(func):
    """Decorator adding a --source option for download source selection."""
    return click.option(
        "--source",
        type=click.Choice(["auto", "osf", "aws"], case_sensitive=False),
        default="auto",
        show_default=True,
        help=(
            "Download source: 'osf' extracts from OSF tar.xz archives, "
            "'aws' fetches individual files from S3, "
            "'auto' estimates which is faster and picks it."
        ),
    )(func)


def format_option(func):
    """Decorator adding a --format option for output format selection."""
    return click.option(
        "--format",
        "output_format",
        type=click.Choice(list(FORMATS), case_sensitive=False),
        default="auto",
        show_default=True,
        help=(
            "Output format: 'table' for human-readable aligned table, "
            "'tsv'/'csv'/'json' for machine-readable output, "
            "'auto' uses table when stdout is a TTY, otherwise tsv."
        ),
    )(func)


def _download_assemblies(
    selected_df: pd.DataFrame,
    file_list_df: pd.DataFrame,
    output: Path,
    cache_dir: Path,
    no_cache: bool,
    threads: int,
    source: str,
) -> list:
    """Choose download source and fetch assemblies."""
    sample_ids = selected_df["sample"].tolist()
    n_genomes = len(sample_ids)

    # Resolve source
    if source == "auto":
        # Figure out how many tarballs we'd need
        tarballs = resolve_tarballs(sample_ids, file_list_df)
        n_tarballs = len(tarballs)
        method, aws_est, osf_est = estimate_download_time(n_genomes, n_tarballs)
        click.echo(
            f"  Estimated time — AWS: ~{aws_est:.0f}s, "
            f"OSF tarballs: ~{osf_est:.0f}s "
            f"({n_tarballs} tarballs)"
        )
        click.echo(f"  Auto-selected source: {method}")
    else:
        method = source

    if method == "aws":
        click.echo("Downloading assemblies from AWS S3...")
        return fetch_from_aws(sample_ids, output, max_workers=min(8, threads + 1))
    else:
        click.echo("Downloading assemblies via OSF tarballs...")
        return fetch_assemblies(
            selected_df, file_list_df, output, cache_dir, no_cache, threads=threads
        )


def _suggest_species(query: str, species_list: list[str]) -> list[str]:
    """Return close species name matches using difflib.

    Parameters
    ----------
    query : str
        The species name the user provided.
    species_list : list[str]
        Available species names to match against.

    Returns
    -------
    list[str]
        Up to 5 close matches.
    """
    return difflib.get_close_matches(query, species_list, n=5, cutoff=0.5)


def _load_toml_filters(filter_file: Path | None) -> dict:
    """Load filter criteria from a TOML file.

    Parameters
    ----------
    filter_file : Path or None
        Path to a TOML filter file.  Returns empty dict if None or missing.

    Returns
    -------
    dict
        Flattened filter dict (keys match CLI option names without dashes).
    """
    if filter_file is None or not filter_file.exists():
        return {}
    with open(filter_file, "rb") as f:
        data = tomllib.load(f)
    # Flatten one level: {"filters": {"species": "..."}} -> {"species": "..."}
    flat: dict = {}
    for section_vals in data.values():
        if isinstance(section_vals, dict):
            flat.update(section_vals)
        # ignore non-dict top-level values
    return flat


# -- Main CLI group --


@click.group()
@click.version_option()
def main():
    """atbfetcher — Fetch genomes from AllTheBacteria for benchmark datasets."""


# -- species subcommand --


@main.command()
@click.argument("species_name")
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for downloaded assemblies.",
)
@click.option("--n", "-n", default=1000, show_default=True, help="Number of genomes to select.")
@click.option("--seed", default=42, show_default=True, help="Random seed for reproducibility.")
@threads_option
@source_option
@quality_filter_option
@cache_options
@verbose_option
def species(
    species_name,
    output,
    n,
    seed,
    threads,
    source,
    quality_filter,
    cache_dir,
    no_cache,
    refresh,
    verbose,
):
    """Fetch a stratified subsample of genomes for a species.

    SPECIES_NAME is the species to fetch (e.g. "Escherichia coli").
    """
    _setup_logging(verbose)

    use_hq = quality_filter != "none"
    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)

    click.echo(f"Loading metadata for {species_name}...")
    species_calls_df = cache.load_species_calls(hq_only=use_hq)
    checkm2_df = cache.load_checkm2()
    file_list_df = cache.load_file_list()

    click.echo(f"Finding samples for {species_name}...")
    samples_df = get_samples_for_species(species_name, species_calls_df)
    click.echo(f"  Found {len(samples_df)} samples")

    if samples_df.empty:
        click.echo(f"No samples found for species: {species_name}", err=True)
        # Feature 1: species typo suggestions
        available = list_species(species_calls_df)
        suggestions = _suggest_species(species_name, available)
        if suggestions:
            click.echo("Did you mean one of:", err=True)
            for s in suggestions:
                click.echo(f"  {s}", err=True)
        sys.exit(1)

    if quality_filter == "qualibact":
        click.echo("Loading Qualibact cutoffs...")
        qualibact_cutoffs = load_qualibact_cutoffs()
        click.echo("Filtering by quality (Qualibact)...")
        hq_df = filter_by_quality(samples_df, checkm2_df, species_name, qualibact_cutoffs)
        click.echo(f"  {len(hq_df)} samples after Qualibact filtering")
    elif quality_filter == "atb":
        click.echo("Using ATB HQ filter (pre-applied)...")
        hq_df = samples_df.merge(checkm2_df, on="sample", how="inner")
        click.echo(f"  {len(hq_df)} high-quality samples")
    else:
        click.echo("Quality filtering: disabled")
        hq_df = samples_df.merge(checkm2_df, on="sample", how="inner")
        click.echo(f"  {len(hq_df)} samples (no quality filter)")

    if hq_df.empty:
        click.echo("No samples passed quality filters.", err=True)
        sys.exit(1)

    click.echo(f"Selecting {n} genomes via stratified sampling...")
    selected_df = stratified_sample(hq_df, n=n, seed=seed)
    click.echo(f"  Selected {len(selected_df)} genomes")

    output.mkdir(parents=True, exist_ok=True)
    _write_summary(selected_df, output, species_name)

    click.echo("Generating selection plot...")
    plot_selection(hq_df, selected_df, output, species_name)

    extracted = _download_assemblies(
        selected_df, file_list_df, output, cache_dir, no_cache, threads, source
    )
    click.echo(f"Done! {len(extracted)} assemblies saved to {output}")


# -- mlst subcommand --


@main.command()
@click.argument("species_name")
@click.option("--scheme", default=None, help="MLST scheme (auto-detected if not specified).")
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for downloaded assemblies.",
)
@click.option("--n", "-n", default=1000, show_default=True, help="Number of genomes to select.")
@click.option("--seed", default=42, show_default=True, help="Random seed for reproducibility.")
@click.option(
    "--strategy",
    type=click.Choice(list(STRATEGIES), case_sensitive=False),
    default="frequency",
    show_default=True,
    help="ST selection strategy.",
)
@threads_option
@source_option
@quality_filter_option
@cache_options
@verbose_option
def mlst(
    species_name,
    scheme,
    output,
    n,
    seed,
    strategy,
    threads,
    source,
    quality_filter,
    cache_dir,
    no_cache,
    refresh,
    verbose,
):
    """Fetch genomes selected by MLST sequence types.

    SPECIES_NAME is the species to fetch (e.g. "Escherichia coli").
    """
    _setup_logging(verbose)

    use_hq = quality_filter != "none"
    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)

    click.echo(f"Loading metadata for {species_name}...")
    species_calls_df = cache.load_species_calls(hq_only=use_hq)
    checkm2_df = cache.load_checkm2()
    file_list_df = cache.load_file_list()

    click.echo(f"Finding samples for {species_name}...")
    samples_df = get_samples_for_species(species_name, species_calls_df)
    click.echo(f"  Found {len(samples_df)} samples")

    if samples_df.empty:
        click.echo(f"No samples found for species: {species_name}", err=True)
        # Feature 1: species typo suggestions
        available = list_species(species_calls_df)
        suggestions = _suggest_species(species_name, available)
        if suggestions:
            click.echo("Did you mean one of:", err=True)
            for s in suggestions:
                click.echo(f"  {s}", err=True)
        sys.exit(1)

    if quality_filter == "qualibact":
        click.echo("Loading Qualibact cutoffs...")
        qualibact_cutoffs = load_qualibact_cutoffs()
        click.echo("Filtering by quality (Qualibact)...")
        hq_df = filter_by_quality(samples_df, checkm2_df, species_name, qualibact_cutoffs)
        click.echo(f"  {len(hq_df)} samples after Qualibact filtering")
    elif quality_filter == "atb":
        click.echo("Using ATB HQ filter (pre-applied)...")
        hq_df = samples_df.merge(checkm2_df, on="sample", how="inner")
        click.echo(f"  {len(hq_df)} high-quality samples")
    else:
        click.echo("Quality filtering: disabled")
        hq_df = samples_df.merge(checkm2_df, on="sample", how="inner")
        click.echo(f"  {len(hq_df)} samples (no quality filter)")

    if hq_df.empty:
        click.echo("No samples passed quality filters.", err=True)
        sys.exit(1)

    click.echo("Loading MLST data...")
    mlst_df = cache.load_mlst()
    suspect_df = load_suspect_contaminations()

    scheme_label = scheme or "auto-detect"
    click.echo(f"Selecting {n} genomes by MLST (scheme: {scheme_label})...")
    selected_df = filter_by_mlst(
        hq_df, mlst_df, scheme=scheme, n=n, suspect_df=suspect_df, seed=seed, strategy=strategy
    )
    click.echo(f"  Selected {len(selected_df)} genomes")

    if selected_df.empty:
        click.echo("No genomes selected by MLST criteria.", err=True)
        sys.exit(1)

    output.mkdir(parents=True, exist_ok=True)
    _write_summary(selected_df, output, species_name)

    click.echo("Generating selection plot...")
    plot_selection(hq_df, selected_df, output, species_name)

    extracted = _download_assemblies(
        selected_df, file_list_df, output, cache_dir, no_cache, threads, source
    )
    click.echo(f"Done! {len(extracted)} assemblies saved to {output}")


# -- accessions subcommand --


@main.command()
@click.argument("accessions_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory for downloaded assemblies.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print what would be downloaded without fetching.",
)
@threads_option
@source_option
@cache_options
@verbose_option
def accessions(
    accessions_file, output, dry_run, threads, source, cache_dir, no_cache, refresh, verbose
):
    """Fetch assemblies for a list of accession IDs.

    ACCESSIONS_FILE is a text file with one accession per line.
    """
    _setup_logging(verbose)

    # Read accession list
    acc_list = [
        line.strip()
        for line in accessions_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    click.echo(f"Read {len(acc_list)} accessions from {accessions_file}")

    # Feature 4: dry-run
    if dry_run:
        click.echo(f"[dry-run] Would download {len(acc_list)} assemblies to {output}:")
        for acc in acc_list:
            click.echo(f"  {acc}")
        return

    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)

    click.echo("Loading file list...")
    file_list_df = cache.load_file_list()

    selected_df = pd.DataFrame({"sample": acc_list})

    extracted = _download_assemblies(
        selected_df, file_list_df, output, cache_dir, no_cache, threads, source
    )
    click.echo(f"Done! {len(extracted)} assemblies saved to {output}")


# -- list-species subcommand --


@main.command("list-species")
@click.option(
    "--raw", is_flag=True, default=False, help="Print original GTDB names without cleaning."
)
@click.option(
    "--count",
    is_flag=True,
    default=False,
    help="Show genome count per species, ordered by count (highest first).",
)
@format_option
@cache_options
@verbose_option
def list_species_cmd(raw, count, output_format, cache_dir, no_cache, refresh, verbose):
    """List all available species in AllTheBacteria."""
    _setup_logging(verbose)

    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)
    species_calls_df = cache.load_species_calls()

    if count:
        # Exclude placeholders then count, ordered highest to lowest
        named = species_calls_df[
            species_calls_df["species"].apply(
                lambda x: pd.notna(x) and not is_placeholder_species(str(x))
            )
        ]
        if not raw:
            named = named.copy()
            named["species"] = named["species"].apply(lambda x: clean_species_name(str(x)))
        counts = named["species"].value_counts().reset_index()
        counts.columns = ["species", "count"]
        fmt = resolve_format(output_format)
        if fmt == "table":
            for _, row in counts.iterrows():
                click.echo(f"{row['count']:>10,}  {row['species']}")
        else:
            print_dataframe(counts, output_format)
    else:
        names = list_species(species_calls_df, raw=raw)
        fmt = resolve_format(output_format)
        if fmt in ("tsv", "csv", "json"):
            df = pd.DataFrame({"species": names})
            print_dataframe(df, output_format)
        else:
            for name in names:
                click.echo(name)


# -- species-count subcommand --


@main.command("species-count")
@click.option(
    "--top", default=0, show_default=True, help="Show only the top N species by count (0 = all)."
)
@format_option
@cache_options
@verbose_option
def species_count(top, output_format, cache_dir, no_cache, refresh, verbose):
    """Show the number of HQ genomes per species in AllTheBacteria."""
    _setup_logging(verbose)

    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)
    species_calls_df = cache.load_species_calls()

    # Exclude GTDB placeholder species (e.g. "Genus sp000746275")
    named = species_calls_df[
        species_calls_df["species"].apply(
            lambda x: pd.notna(x) and not is_placeholder_species(str(x))
        )
    ]
    counts = named["species"].value_counts().reset_index()
    counts.columns = ["species", "count"]

    if top > 0:
        counts = counts.head(top)

    fmt = resolve_format(output_format)
    if fmt in ("tsv", "csv", "json"):
        print_dataframe(counts, output_format)
    else:
        click.echo(f"{'Species':<50} {'Count':>10}")
        click.echo("-" * 62)
        for _, row in counts.iterrows():
            click.echo(f"{row['species']:<50} {row['count']:>10,}")
        click.echo(f"\nTotal: {counts['count'].sum():,} HQ genomes across {len(counts)} species")


# -- query subcommand --


@main.command()
@click.option("--species", "-s", default=None, help="Filter by species (e.g. 'Escherichia coli').")
@click.option(
    "--country", "-c", default=None, help="Filter by country (prefix match, e.g. 'United Kingdom')."
)
@click.option("--year-from", type=int, default=None, help="Minimum collection year (inclusive).")
@click.option("--year-to", type=int, default=None, help="Maximum collection year (inclusive).")
@click.option("--host", default=None, help="Filter by host organism (e.g. 'Homo sapiens').")
@click.option(
    "--isolation-source", default=None, help="Filter by isolation source (e.g. 'blood', 'stool')."
)
@click.option(
    "--hq-only/--no-hq",
    default=True,
    show_default=True,
    help="Only include ATB high-quality genomes.",
)
@click.option(
    "--min-completeness", type=float, default=None, help="Minimum CheckM2 completeness (e.g. 95.0)."
)
@click.option(
    "--max-contamination",
    type=float,
    default=None,
    help="Maximum CheckM2 contamination (e.g. 5.0).",
)
@click.option("--min-genome-size", type=int, default=None, help="Minimum genome size in bases.")
@click.option("--max-genome-size", type=int, default=None, help="Maximum genome size in bases.")
@click.option(
    "--n",
    "-n",
    type=int,
    default=None,
    help="Maximum number of genomes to return (random subsample).",
)
@click.option(
    "--seed", default=42, show_default=True, help="Random seed for reproducible subsampling."
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory for downloaded assemblies. If not set, only lists accessions.",
)
@click.option(
    "--db-path",
    type=click.Path(path_type=Path),
    default=None,
    help=f"Path to ATB SQLite database. Default: looks in cache dir for {SQLITE_FILENAME}.",
)
@click.option(
    "--filter",
    "filter_file",
    type=click.Path(path_type=Path),
    default=None,
    help="TOML file with filter criteria (CLI flags override TOML values).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print what would be downloaded without fetching.",
)
@format_option
@threads_option
@source_option
@cache_options
@verbose_option
def query(
    species,
    country,
    year_from,
    year_to,
    host,
    isolation_source,
    hq_only,
    min_completeness,
    max_contamination,
    min_genome_size,
    max_genome_size,
    n,
    seed,
    output,
    db_path,
    filter_file,
    dry_run,
    output_format,
    threads,
    source,
    cache_dir,
    no_cache,
    refresh,
    verbose,
):
    """Query the ATB metadata database to select genomes by metadata.

    Uses the ATB SQLite metadata database to filter genomes by species,
    country, collection date, host, isolation source, and quality metrics.

    \b
    Examples:
      # List E. coli from the UK collected 2020-2023
      atbfetcher query --species "Escherichia coli" --country "United Kingdom" \\
        --year-from 2020 --year-to 2023

    \b
      # Fetch 500 S. aureus from blood samples
      atbfetcher query --species "Staphylococcus aureus" \\
        --isolation-source blood --n 500 --output ./saureus_blood

    \b
      # Use a TOML filter file
      atbfetcher query --filter my_query.toml

    \b
    The SQLite database must be downloaded first:
      atbfetcher download-db
    """
    _setup_logging(verbose)

    # Feature 6: load TOML filter file, CLI flags take precedence (override)
    toml_filters = _load_toml_filters(filter_file)
    if toml_filters:
        click.echo(f"Loaded filter file: {filter_file}", err=True)

    def _resolve(cli_val, toml_key, default=None):
        """Return CLI value if explicitly set, else fall back to TOML, then default."""
        if cli_val is not None:
            return cli_val
        return toml_filters.get(toml_key, default)

    # Apply TOML defaults where CLI didn't provide values
    species = _resolve(species, "species")
    country = _resolve(country, "country")
    year_from = _resolve(year_from, "year_from")
    year_to = _resolve(year_to, "year_to")
    host = _resolve(host, "host")
    isolation_source = _resolve(isolation_source, "isolation_source")
    if "hq_only" in toml_filters and hq_only is True:
        hq_only = toml_filters["hq_only"]
    min_completeness = _resolve(min_completeness, "min_completeness")
    max_contamination = _resolve(max_contamination, "max_contamination")
    min_genome_size = _resolve(min_genome_size, "min_genome_size")
    max_genome_size = _resolve(max_genome_size, "max_genome_size")
    n = _resolve(n, "n")

    # Locate the SQLite database
    db = find_sqlite_db(db_path, Path(cache_dir))
    if db is None:
        click.echo(
            "ATB SQLite database not found.\n\n"
            "Download it with:\n"
            "  atbfetcher download-db\n\n"
            "Or specify its location with --db-path.",
            err=True,
        )
        sys.exit(1)

    click.echo(f"Using database: {db}")

    # Run query
    click.echo("Querying metadata...")
    results = query_metadata(
        db,
        species=species,
        country=country,
        year_from=year_from,
        year_to=year_to,
        host=host,
        isolation_source=isolation_source,
        hq_only=hq_only,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        min_genome_size=min_genome_size,
        max_genome_size=max_genome_size,
        limit=n,
        seed=seed,
    )

    if results.empty:
        click.echo("No genomes matched the query filters.", err=True)
        # Feature 1: species suggestions when species filter was set
        if species:
            try:
                conn = sqlite3.connect(str(db))
                rows = conn.execute(
                    "SELECT DISTINCT sylph_species FROM assembly "
                    "WHERE sylph_species IS NOT NULL"
                ).fetchall()
                conn.close()
                all_species = [r[0] for r in rows]
                suggestions = _suggest_species(species, all_species)
                if suggestions:
                    click.echo("Did you mean one of:", err=True)
                    for s in suggestions:
                        click.echo(f"  {s}", err=True)
            except Exception:
                pass
        sys.exit(1)

    click.echo(f"  Found {len(results)} matching genomes")

    # Print summary of filters applied
    filters = []
    if species:
        filters.append(f"species={species}")
    if country:
        filters.append(f"country={country}")
    if year_from or year_to:
        yr = f"{year_from or '...'}-{year_to or '...'}"
        filters.append(f"years={yr}")
    if host:
        filters.append(f"host={host}")
    if isolation_source:
        filters.append(f"isolation_source={isolation_source}")
    if hq_only:
        filters.append("hq_only=True")
    if filters:
        click.echo(f"  Filters: {', '.join(filters)}")

    if output:
        # Feature 4: dry-run for downloads
        if dry_run:
            sample_ids = results["sample"].tolist()
            click.echo(
                f"[dry-run] Would download {len(sample_ids)} assemblies to {output}:"
            )
            for sid in sample_ids:
                click.echo(f"  {sid}")
            return

        # Download assemblies
        output.mkdir(parents=True, exist_ok=True)

        # Save query results as summary TSV
        safe_name = (species or "query").replace(" ", "_")
        summary_path = output / f"{safe_name}_query_results.tsv"
        results.to_csv(summary_path, sep="\t", index=False)
        click.echo(f"  Query results saved to {summary_path}")

        # Download using AWS (query results have aws_url directly)
        sample_ids = results["sample"].tolist()
        n_genomes = len(sample_ids)

        if source == "auto":
            # For query mode, prefer AWS since we don't have tarball info
            # readily available from the SQLite DB
            cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)
            try:
                file_list_df = cache.load_file_list()
                tarballs = resolve_tarballs(sample_ids, file_list_df)
                n_tarballs = len(tarballs)
                method, aws_est, osf_est = estimate_download_time(n_genomes, n_tarballs)
                click.echo(
                    f"  Estimated time — AWS: ~{aws_est:.0f}s, "
                    f"OSF tarballs: ~{osf_est:.0f}s ({n_tarballs} tarballs)"
                )
                click.echo(f"  Auto-selected source: {method}")
            except Exception:
                method = "aws"
                file_list_df = None
                click.echo("  Using AWS (file list unavailable for tarball estimation)")
        else:
            method = source
            file_list_df = None

        if method == "aws":
            click.echo("Downloading assemblies from AWS S3...")
            extracted = fetch_from_aws(sample_ids, output, max_workers=min(8, threads + 1))
        else:
            if file_list_df is None:
                cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)
                file_list_df = cache.load_file_list()
            click.echo("Downloading assemblies via OSF tarballs...")
            selected_df = pd.DataFrame({"sample": sample_ids})
            extracted = fetch_assemblies(
                selected_df, file_list_df, output, cache_dir, no_cache, threads=threads
            )

        click.echo(f"Done! {len(extracted)} assemblies saved to {output}")
    else:
        # List mode: print results using selected output format
        print_dataframe(results[["sample"]], output_format)


# -- query list-countries subcommand --


@main.command("list-countries")
@click.option(
    "--species", "-s", default=None, help="Only show countries with genomes of this species."
)
@click.option(
    "--db-path", type=click.Path(path_type=Path), default=None, help="Path to ATB SQLite database."
)
@cache_options
@verbose_option
def list_countries_cmd(species, db_path, cache_dir, no_cache, refresh, verbose):
    """List available countries in the ATB metadata database."""
    _setup_logging(verbose)

    db = find_sqlite_db(db_path, Path(cache_dir))
    if db is None:
        click.echo(
            "ATB SQLite database not found. Run 'atbfetcher download-db' first.",
            err=True,
        )
        sys.exit(1)

    countries = list_countries(db, species=species)
    for c in countries:
        click.echo(c)
    click.echo(f"\n{len(countries)} countries")


# -- query list-hosts subcommand --


@main.command("list-hosts")
@click.option("--species", "-s", default=None, help="Only show hosts with genomes of this species.")
@click.option(
    "--db-path", type=click.Path(path_type=Path), default=None, help="Path to ATB SQLite database."
)
@cache_options
@verbose_option
def list_hosts_cmd(species, db_path, cache_dir, no_cache, refresh, verbose):
    """List available host organisms in the ATB metadata database."""
    _setup_logging(verbose)

    db = find_sqlite_db(db_path, Path(cache_dir))
    if db is None:
        click.echo(
            "ATB SQLite database not found. Run 'atbfetcher download-db' first.",
            err=True,
        )
        sys.exit(1)

    hosts = list_hosts(db, species=species)
    for h in hosts:
        click.echo(h)
    click.echo(f"\n{len(hosts)} hosts")


# -- download-db subcommand --


@main.command("download-db")
@cache_options
@verbose_option
def download_db_cmd(cache_dir, no_cache, refresh, verbose):
    """Download the ATB SQLite metadata database.

    Downloads the compressed database (~2 GB) from OSF and decompresses it
    (~27 GB). Required for the ``query``, ``list-countries``, and
    ``list-hosts`` commands.

    The database is stored in the cache directory (default: ~/.atbfetcher/).
    """
    _setup_logging(verbose)

    cache_dir = Path(cache_dir)
    db_path = cache_dir / SQLITE_FILENAME

    if db_path.exists() and not refresh:
        click.echo(f"Database already exists at {db_path}")
        click.echo(f"  Size: {db_path.stat().st_size / 1e9:.1f} GB")
        click.echo("Use --refresh to re-download.")
        return

    if db_path.exists() and refresh:
        click.echo("Removing existing database for re-download...")
        db_path.unlink()

    click.echo("Downloading ATB SQLite metadata database...")
    click.echo(f"  Source: {SQLITE_URL}")
    click.echo(f"  Destination: {cache_dir}/")
    click.echo("  Download size: ~2 GB, uncompressed: ~27 GB")
    click.echo()

    try:
        result_path = download_sqlite_db(cache_dir)
        click.echo()
        click.echo(f"Database ready at {result_path}")
        click.echo(f"  Size: {result_path.stat().st_size / 1e9:.1f} GB")
    except RuntimeError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


# -- info subcommand (Feature 3) --


@main.command()
@click.argument("sample_accession")
@click.option(
    "--db-path",
    type=click.Path(path_type=Path),
    default=None,
    help=f"Path to ATB SQLite database. Default: looks in cache dir for {SQLITE_FILENAME}.",
)
@cache_options
@verbose_option
def info(sample_accession, db_path, cache_dir, no_cache, refresh, verbose):
    """Show a detailed summary for a single sample accession.

    \b
    Example:
      atbfetcher info SAMN12345678
    """
    _setup_logging(verbose)

    db = find_sqlite_db(db_path, Path(cache_dir))
    if db is None:
        click.echo(
            "ATB SQLite database not found. Run 'atbfetcher download-db' first.",
            err=True,
        )
        sys.exit(1)

    conn = sqlite3.connect(str(db))
    try:
        # Assembly info
        asm_row = conn.execute(
            "SELECT sample_accession, sylph_species, hq_filter, aws_url, "
            "asm_fasta_on_osf FROM assembly WHERE sample_accession = ?",
            [sample_accession],
        ).fetchone()

        if asm_row is None:
            click.echo(f"Sample not found: {sample_accession}", err=True)
            sys.exit(1)

        click.echo(f"\n{'='*60}")
        click.echo(f"  Sample: {asm_row[0]}")
        click.echo(f"{'='*60}")

        click.echo("\n[Assembly]")
        click.echo(f"  Species      : {asm_row[1] or 'N/A'}")
        click.echo(f"  HQ filter    : {asm_row[2] or 'N/A'}")
        click.echo(f"  On OSF       : {'Yes' if asm_row[4] else 'No'}")
        click.echo(f"  AWS URL      : {asm_row[3] or 'N/A'}")

        # CheckM2 stats
        chk_row = conn.execute(
            "SELECT Completeness_Specific, Contamination, Genome_Size, "
            "GC_Content, Contig_N50 FROM checkm2 WHERE sample_accession = ?",
            [sample_accession],
        ).fetchone()

        click.echo("\n[Assembly Statistics (CheckM2)]")
        if chk_row:
            click.echo(f"  Completeness : {chk_row[0]:.1f}%")
            click.echo(f"  Contamination: {chk_row[1]:.2f}%")
            click.echo(f"  Genome size  : {chk_row[2]:,} bp" if chk_row[2] else "  Genome size  : N/A")
            click.echo(f"  GC content   : {chk_row[3]:.1f}%" if chk_row[3] else "  GC content   : N/A")
            click.echo(f"  Contig N50   : {chk_row[4]:,} bp" if chk_row[4] else "  Contig N50   : N/A")
        else:
            click.echo("  (no CheckM2 data)")

        # MLST — may have multiple scheme rows; table may not exist in all DBs
        mlst_rows = []
        try:
            mlst_rows = conn.execute(
                "SELECT mlst_scheme, mlst_st, mlst_status FROM mlst "
                "WHERE sample_accession = ? ORDER BY mlst_scheme",
                [sample_accession],
            ).fetchall()
        except Exception:
            # Try alternative column names used in some ATB schema versions
            try:
                mlst_rows = conn.execute(
                    "SELECT scheme, st, status FROM mlst "
                    "WHERE sample = ? ORDER BY scheme",
                    [sample_accession],
                ).fetchall()
            except Exception:
                mlst_rows = []

        click.echo("\n[MLST]")
        if mlst_rows:
            for row in mlst_rows:
                click.echo(f"  Scheme: {row[0]}  ST: {row[1]}  Status: {row[2]}")
        else:
            click.echo("  (no MLST data in database)")

        # ENA metadata
        ena_rows = conn.execute(
            "SELECT e.country, e.collection_date, e.host, e.isolation_source, "
            "r.run_accession "
            "FROM run r JOIN ena_202505_used e ON r.run_accession = e.run_accession "
            "WHERE r.sample_accession = ?",
            [sample_accession],
        ).fetchall()

        click.echo("\n[ENA Metadata]")
        if ena_rows:
            for row in ena_rows:
                click.echo(f"  Run              : {row[4]}")
                click.echo(f"  Country          : {row[0] or 'N/A'}")
                click.echo(f"  Collection date  : {row[1] or 'N/A'}")
                click.echo(f"  Host             : {row[2] or 'N/A'}")
                click.echo(f"  Isolation source : {row[3] or 'N/A'}")
                if len(ena_rows) > 1:
                    click.echo()
        else:
            click.echo("  (no ENA metadata)")

        click.echo()
    finally:
        conn.close()


# -- summarise subcommand (Feature 7) --


@main.command()
@click.option(
    "--by",
    "-b",
    default="species",
    show_default=True,
    help="Column to group by (e.g. species, country, host, isolation_source).",
)
@click.option(
    "--top",
    "-t",
    default=0,
    show_default=True,
    help="Show only the top N groups by count (0 = all).",
)
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(path_type=Path),
    default=None,
    help="Input TSV file (e.g. from query output). Reads from stdin if not set.",
)
@click.option(
    "--db-path",
    type=click.Path(path_type=Path),
    default=None,
    help=f"Path to ATB SQLite database for direct querying.",
)
@click.option("--species", "-s", default=None, help="Pre-filter by species before summarising.")
@format_option
@cache_options
@verbose_option
def summarise(by, top, input_file, db_path, species, output_format, cache_dir, no_cache, refresh, verbose):
    """Compute group-by counts and stats from ATB metadata.

    Can read from a TSV file, stdin, or query the SQLite database directly.

    \b
    Examples:
      # Top 20 species in the database
      atbfetcher summarise --by species --top 20

      # Summarise piped query results by country
      atbfetcher query --species "Escherichia coli" | atbfetcher summarise --by country --input -

      # Summarise from a saved TSV
      atbfetcher summarise --by isolation_source --input results.tsv
    """
    _setup_logging(verbose)

    df: pd.DataFrame | None = None

    # Determine data source
    if input_file is not None:
        if str(input_file) == "-":
            # Read from stdin
            import io as _io
            raw = sys.stdin.read()
            if raw.strip():
                df = pd.read_csv(_io.StringIO(raw), sep="\t")
        else:
            df = pd.read_csv(input_file, sep="\t")
    elif not sys.stdin.isatty():
        # Piped input — only read if stdin has content
        import io as _io
        raw = sys.stdin.read()
        if raw.strip():
            df = pd.read_csv(_io.StringIO(raw), sep="\t")

    if df is None:
        # Query the SQLite database
        db = find_sqlite_db(db_path, Path(cache_dir))
        if db is None:
            click.echo(
                "No input provided and ATB SQLite database not found.\n\n"
                "Either pipe query output, provide --input, or run 'atbfetcher download-db'.",
                err=True,
            )
            sys.exit(1)

        click.echo(f"Using database: {db}", err=True)
        df = query_metadata(
            db,
            species=species,
            country=None,
            hq_only=True,
        )

    if df.empty:
        click.echo("No data to summarise.", err=True)
        sys.exit(1)

    if by not in df.columns:
        available_cols = ", ".join(df.columns.tolist())
        click.echo(
            f"Column '{by}' not found in data. Available columns: {available_cols}",
            err=True,
        )
        sys.exit(1)

    # Group by the chosen column and compute counts + numeric stats
    grouped = df.groupby(by)
    count_series = grouped.size().rename("count").reset_index()
    count_series = count_series.sort_values("count", ascending=False).reset_index(drop=True)

    # Add mean stats for numeric columns if present
    numeric_cols = [
        c for c in ["Completeness_Specific", "Contamination", "Genome_Size", "Contig_N50"]
        if c in df.columns
    ]
    if numeric_cols:
        agg = grouped[numeric_cols].mean().round(2).reset_index()
        count_series = count_series.merge(agg, on=by, how="left")

    if top > 0:
        count_series = count_series.head(top)

    print_dataframe(count_series, output_format)


# -- mlst-query subcommand (Feature 5) --


@main.command("mlst-query")
@click.option(
    "--species", "-s", default=None, help="Filter by species (e.g. 'Escherichia coli')."
)
@click.option("--scheme", default=None, help="Filter by MLST scheme name (e.g. 'ecoli_achtman_4').")
@click.option("--st", default=None, help="Filter by sequence type number (e.g. '131').")
@click.option(
    "--hq-only/--no-hq",
    default=True,
    show_default=True,
    help="Only include ATB high-quality genomes.",
)
@click.option(
    "--n",
    "-n",
    type=int,
    default=None,
    help="Maximum number of results to return.",
)
@format_option
@cache_options
@verbose_option
def mlst_query(species, scheme, st, hq_only, n, output_format, cache_dir, no_cache, refresh, verbose):
    """Query samples by MLST data — output accessions + MLST info.

    Filters the MLST dataset by ST number, scheme name, and/or species.
    Prints a table of matching sample accessions with their MLST data.

    \b
    Examples:
      # All E. coli ST131 samples
      atbfetcher mlst-query --species "Escherichia coli" --st 131

      # Samples typed with the ecoli_achtman_4 scheme
      atbfetcher mlst-query --scheme ecoli_achtman_4

      # E. coli ST131 output as TSV
      atbfetcher mlst-query --species "Escherichia coli" --st 131 --format tsv
    """
    _setup_logging(verbose)

    cache = MetadataCache(cache_dir=cache_dir, no_cache=no_cache, refresh=refresh)

    click.echo("Loading MLST data...", err=True)
    mlst_df = cache.load_mlst()

    # Filter by scheme
    if scheme:
        mlst_df = mlst_df[mlst_df["mlst_scheme"] == scheme]
        if mlst_df.empty:
            click.echo(f"No results for scheme '{scheme}'.", err=True)
            sys.exit(1)

    # Filter by ST
    if st:
        mlst_df = mlst_df[mlst_df["mlst_st"].astype(str) == str(st)]
        if mlst_df.empty:
            click.echo(f"No results for ST '{st}'.", err=True)
            sys.exit(1)

    # Filter by species (requires species calls)
    if species:
        click.echo("Loading species calls...", err=True)
        species_calls_df = cache.load_species_calls(hq_only=hq_only)
        samples_df = get_samples_for_species(species, species_calls_df)
        if samples_df.empty:
            click.echo(f"No samples found for species: {species}", err=True)
            available = list_species(species_calls_df)
            suggestions = _suggest_species(species, available)
            if suggestions:
                click.echo("Did you mean one of:", err=True)
                for s in suggestions:
                    click.echo(f"  {s}", err=True)
            sys.exit(1)
        mlst_df = mlst_df[mlst_df["sample"].isin(samples_df["sample"])]
        if mlst_df.empty:
            click.echo(f"No MLST results for species '{species}'.", err=True)
            sys.exit(1)
    elif hq_only:
        # Apply HQ filter without species constraint
        click.echo("Loading species calls for HQ filter...", err=True)
        species_calls_df = cache.load_species_calls(hq_only=True)
        mlst_df = mlst_df[mlst_df["sample"].isin(species_calls_df["sample"])]

    if n is not None:
        mlst_df = mlst_df.head(n)

    click.echo(f"Found {len(mlst_df)} matching samples.", err=True)
    print_dataframe(mlst_df, output_format)


# -- config subcommand (Feature 8) --


@main.group()
def config():
    """Manage atbfetcher user configuration.

    Settings in ``~/.atbfetcher/config.toml`` serve as defaults that CLI flags
    can always override.

    \b
    Subcommands:
      init   Create a default config file
      show   Display the current config
      set    Set a config value
    """


@config.command("init")
@click.option("--force", is_flag=True, default=False, help="Overwrite existing config.")
def config_init(force):
    """Create a default config file at ~/.atbfetcher/config.toml."""
    if CONFIG_PATH.exists() and not force:
        click.echo(f"Config already exists at {CONFIG_PATH}")
        click.echo("Use --force to overwrite.")
        return

    cfg = default_config()
    write_config(cfg)
    click.echo(f"Created default config at {CONFIG_PATH}")
    _print_config_dict(cfg)


@config.command("show")
def config_show():
    """Display the current configuration."""
    cfg = load_config()
    if not cfg:
        click.echo(f"No config file found at {CONFIG_PATH}")
        click.echo("Run 'atbfetcher config init' to create one.")
        return

    click.echo(f"Config file: {CONFIG_PATH}\n")
    _print_config_dict(cfg)


@config.command("set")
@click.argument("key")
@click.argument("value")
@click.option(
    "--section",
    default="defaults",
    show_default=True,
    help="Config section to write to.",
)
def config_set(key, value, section):
    """Set a config value.

    \b
    Examples:
      atbfetcher config set cache_dir /data/atb_cache
      atbfetcher config set output_format tsv
      atbfetcher config set threads 8
    """
    cfg = load_config()
    if not cfg:
        cfg = default_config()

    if section not in cfg:
        cfg[section] = {}

    # Try to coerce numeric values
    coerced: str | int | float = value
    try:
        coerced = int(value)
    except ValueError:
        try:
            coerced = float(value)
        except ValueError:
            coerced = value

    cfg[section][key] = coerced
    write_config(cfg)
    click.echo(f"Set [{section}] {key} = {coerced!r}")
    click.echo(f"Config saved to {CONFIG_PATH}")


def _print_config_dict(cfg: dict) -> None:
    """Pretty-print a config dict."""
    for section, values in cfg.items():
        click.echo(f"[{section}]")
        if isinstance(values, dict):
            for k, v in values.items():
                click.echo(f"  {k} = {v!r}")
        else:
            click.echo(f"  {values!r}")
        click.echo()


if __name__ == "__main__":
    main()
