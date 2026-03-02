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
from atbfetcher.mlst import filter_by_mlst, load_suspect_contaminations
from atbfetcher.plotting import plot_selection
from atbfetcher.quality import filter_by_quality
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
    "sample", "species", "Completeness_Specific", "Contamination",
    "Genome_Size", "GC_Content", "Contig_N50",
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
        "--verbose", "-v",
        is_flag=True,
        default=False,
        help="Enable verbose (debug) logging.",
    )(func)


def threads_option(func):
    """Decorator adding a --threads option."""
    return click.option(
        "--threads", "-t",
        default=DEFAULT_THREADS,
        show_default=True,
        help="Number of threads for XZ decompression.",
    )(func)


def quality_filter_option(func):
    """Decorator adding a --quality-filter option."""
    return click.option(
        "--quality-filter", "-q",
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
        method, aws_est, osf_est = estimate_download_time(
            n_genomes, n_tarballs
        )
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
@click.option("--output", "-o", required=True, type=click.Path(path_type=Path),
              help="Output directory for downloaded assemblies.")
@click.option("--n", "-n", default=1000, show_default=True,
              help="Number of genomes to select.")
@click.option("--seed", default=42, show_default=True,
              help="Random seed for reproducibility.")
@threads_option
@source_option
@quality_filter_option
@cache_options
@verbose_option
def species(species_name, output, n, seed, threads, source, quality_filter,
            cache_dir, no_cache, refresh, verbose):
    """Fetch a stratified subsample of genomes for a species.

    SPECIES_NAME is the species to fetch (e.g. "Escherichia coli").
    """
    _setup_logging(verbose)
    logger = logging.getLogger(__name__)

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
@click.option("--scheme", default=None,
              help="MLST scheme (auto-detected if not specified).")
@click.option("--output", "-o", required=True, type=click.Path(path_type=Path),
              help="Output directory for downloaded assemblies.")
@click.option("--n", "-n", default=1000, show_default=True,
              help="Number of genomes to select.")
@click.option("--seed", default=42, show_default=True,
              help="Random seed for reproducibility.")
@threads_option
@source_option
@quality_filter_option
@cache_options
@verbose_option
def mlst(species_name, scheme, output, n, seed, threads, source, quality_filter,
         cache_dir, no_cache, refresh, verbose):
    """Fetch genomes selected by MLST sequence types.

    SPECIES_NAME is the species to fetch (e.g. "Escherichia coli").
    """
    _setup_logging(verbose)
    logger = logging.getLogger(__name__)

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
    # Look in bundled data/ first, then fall back to cache directory
    data_dir = Path(__file__).resolve().parent.parent.parent / "data"
    mlst_path = data_dir / "mlst_processed_all_samples.tsv.xz"
    if not mlst_path.exists():
        mlst_path = cache_dir / "mlst_processed_all_samples.tsv.xz"
    if not mlst_path.exists():
        click.echo(
            "MLST data not found. Place mlst_processed_all_samples.tsv.xz "
            f"in {data_dir} or {cache_dir}.",
            err=True,
        )
        sys.exit(1)

    mlst_df = pd.read_csv(mlst_path, sep="\t")
    suspect_df = load_suspect_contaminations()

    scheme_label = scheme or "auto-detect"
    click.echo(f"Selecting {n} genomes by MLST (scheme: {scheme_label})...")
    selected_df = filter_by_mlst(
        hq_df, mlst_df, scheme=scheme, n=n, suspect_df=suspect_df, seed=seed
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
@click.option("--output", "-o", required=True, type=click.Path(path_type=Path),
              help="Output directory for downloaded assemblies.")
@threads_option
@source_option
@cache_options
@verbose_option
def accessions(accessions_file, output, threads, source,
               cache_dir, no_cache, refresh, verbose):
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
@click.option("--raw", is_flag=True, default=False,
              help="Print original GTDB names without cleaning.")
@click.option("--count", is_flag=True, default=False,
              help="Show genome count per species, ordered by count (highest first).")
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
            named["species"] = named["species"].apply(
                lambda x: clean_species_name(str(x))
            )
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
@click.option("--top", default=0, show_default=True,
              help="Show only the top N species by count (0 = all).")
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
    counts = (
        named["species"]
        .value_counts()
        .reset_index()
    )
    counts.columns = ["species", "count"]

    if top > 0:
        counts = counts.head(top)

    click.echo(f"{'Species':<50} {'Count':>10}")
    click.echo("-" * 62)
    for _, row in counts.iterrows():
        click.echo(f"{row['species']:<50} {row['count']:>10,}")

    click.echo(f"\nTotal: {counts['count'].sum():,} HQ genomes across {len(counts)} species")


if __name__ == "__main__":
    main()
