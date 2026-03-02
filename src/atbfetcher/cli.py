"""Command-line interface for atbfetcher.

Provides subcommands to fetch genomes from AllTheBacteria:
- ``species``: Subsample by species with stratified sampling
- ``mlst``: Select by MLST sequence types
- ``accessions``: Fetch specific accessions
- ``list-species``: Print available species names
- ``species-count``: Show genome counts per species
"""

import logging
import sys
from pathlib import Path

import click
import pandas as pd
from rich.logging import RichHandler

from atbfetcher.download import (
    DEFAULT_THREADS,
    estimate_download_time,
    fetch_assemblies,
    fetch_from_aws,
    resolve_tarballs,
)
from atbfetcher.metadata import DEFAULT_CACHE_DIR, MetadataCache, load_qualibact_cutoffs
from atbfetcher.mlst import STRATEGIES, filter_by_mlst, load_suspect_contaminations
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
@threads_option
@source_option
@cache_options
@verbose_option
def accessions(accessions_file, output, threads, source, cache_dir, no_cache, refresh, verbose):
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
@cache_options
@verbose_option
def list_species_cmd(raw, count, cache_dir, no_cache, refresh, verbose):
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
        for _, row in counts.iterrows():
            click.echo(f"{row['count']:>10,}  {row['species']}")
    else:
        names = list_species(species_calls_df, raw=raw)
        for name in names:
            click.echo(name)


# -- species-count subcommand --


@main.command("species-count")
@click.option(
    "--top", default=0, show_default=True, help="Show only the top N species by count (0 = all)."
)
@cache_options
@verbose_option
def species_count(top, cache_dir, no_cache, refresh, verbose):
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
    The SQLite database must be downloaded first:
      atbfetcher download-db
    """
    _setup_logging(verbose)

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
        # List mode: print accessions to stdout
        for _, row in results.iterrows():
            click.echo(row["sample"])


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


if __name__ == "__main__":
    main()
