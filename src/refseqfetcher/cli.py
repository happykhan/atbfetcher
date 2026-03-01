"""Command-line interface for refseqfetcher.

Fetches genomes from NCBI RefSeq using the ``datasets`` CLI tool
(installed via conda/pixi from ncbi-datasets-cli).

Subcommands:
- ``species``: Download a stratified subsample for a species
- ``accessions``: Download specific RefSeq accessions
"""

import json
import logging
import subprocess
import sys
from pathlib import Path

import click
import pandas as pd

from atbfetcher.metadata import load_qualibact_cutoffs
from atbfetcher.plotting import plot_selection
from atbfetcher.quality import filter_by_quality
from atbfetcher.sampling import stratified_sample

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging for the CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def _run_datasets_command(args: list[str]) -> str:
    """Run an NCBI datasets CLI command and return stdout.

    Raises
    ------
    click.ClickException
        If the datasets command fails or is not found.
    """
    try:
        result = subprocess.run(
            ["datasets"] + args,
            capture_output=True,
            text=True,
            check=True,
            timeout=300,
        )
        return result.stdout
    except FileNotFoundError:
        raise click.ClickException(
            "NCBI 'datasets' CLI not found. "
            "Install it with: pixi add ncbi-datasets-cli"
        )
    except subprocess.CalledProcessError as e:
        raise click.ClickException(f"datasets command failed: {e.stderr}")


def _get_refseq_metadata(taxon: str) -> pd.DataFrame:
    """Fetch RefSeq genome metadata for a taxon using NCBI datasets.

    Returns a DataFrame with columns: accession, GC_Content, Genome_Size,
    Completeness, Contamination, contig_count.
    """
    output = _run_datasets_command([
        "summary", "genome", "taxon", taxon,
        "--assembly-source", "refseq",
        "--as-json-lines",
    ])

    records = []
    for line in output.strip().splitlines():
        if not line.strip():
            continue
        try:
            data = json.loads(line)
        except json.JSONDecodeError:
            continue

        # Extract fields from nested NCBI datasets JSON
        accession = data.get("accession", "")
        stats = data.get("assembly_stats", {})
        info = data.get("assembly_info", {})
        checkm = info.get("checkm_info", {})

        records.append({
            "sample": accession,
            "Name": accession,
            "accession": accession,
            "GC_Content": stats.get("gc_percent", None),
            "Genome_Size": stats.get("total_sequence_length", None),
            "Completeness": checkm.get("completeness", 100.0),
            "Contamination": checkm.get("contamination", 0.0),
            "contig_count": stats.get("number_of_contigs", None),
        })

    df = pd.DataFrame(records)
    logger.info("Found %d RefSeq assemblies for %s", len(df), taxon)
    return df


def _download_accessions(accessions: list[str], output_dir: Path) -> list[Path]:
    """Download genome assemblies by accession using NCBI datasets.

    Parameters
    ----------
    accessions : list[str]
        RefSeq accession IDs to download.
    output_dir : Path
        Directory to store downloaded files.

    Returns
    -------
    list[Path]
        Paths to downloaded assembly files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    zip_path = output_dir / "ncbi_dataset.zip"

    # Download using datasets CLI
    acc_str = ",".join(accessions)
    _run_datasets_command([
        "download", "genome", "accession", acc_str,
        "--include", "genome",
        "--filename", str(zip_path),
    ])

    # Unzip
    import zipfile
    extracted = []
    with zipfile.ZipFile(zip_path, "r") as zf:
        for name in zf.namelist():
            if name.endswith((".fna", ".fna.gz", ".fa", ".fa.gz")):
                zf.extract(name, output_dir)
                extracted.append(output_dir / name)

    # Clean up zip
    zip_path.unlink(missing_ok=True)

    return extracted


# -- Main CLI group --

@click.group()
@click.version_option(package_name="atbfetcher")
def main():
    """refseqfetcher — Fetch genomes from NCBI RefSeq for benchmark datasets."""


# -- species subcommand --

@main.command()
@click.argument("species_name")
@click.option("--output", "-o", required=True, type=click.Path(path_type=Path),
              help="Output directory for downloaded assemblies.")
@click.option("--n", "-n", default=1000, show_default=True,
              help="Number of genomes to select.")
@click.option("--seed", default=42, show_default=True,
              help="Random seed for reproducibility.")
@click.option("--verbose", "-v", is_flag=True, default=False,
              help="Enable verbose logging.")
def species(species_name, output, n, seed, verbose):
    """Fetch a stratified subsample of RefSeq genomes for a species.

    SPECIES_NAME is the species to fetch (e.g. "Escherichia coli").
    """
    _setup_logging(verbose)

    click.echo(f"Fetching RefSeq metadata for {species_name}...")
    metadata_df = _get_refseq_metadata(species_name)

    if metadata_df.empty:
        click.echo(f"No RefSeq assemblies found for: {species_name}", err=True)
        sys.exit(1)

    click.echo(f"  Found {len(metadata_df)} RefSeq assemblies")

    # Apply quality filtering
    click.echo("Loading Qualibact cutoffs...")
    qualibact_cutoffs = load_qualibact_cutoffs()

    click.echo("Filtering by quality...")
    hq_df = filter_by_quality(metadata_df, metadata_df, species_name, qualibact_cutoffs)
    click.echo(f"  {len(hq_df)} high-quality assemblies")

    if hq_df.empty:
        click.echo("No assemblies passed quality filters.", err=True)
        sys.exit(1)

    click.echo(f"Selecting {n} genomes via stratified sampling...")
    selected_df = stratified_sample(hq_df, n=n, seed=seed)
    click.echo(f"  Selected {len(selected_df)} genomes")

    click.echo("Generating selection plot...")
    plot_selection(hq_df, selected_df, Path(output), species_name)

    click.echo("Downloading assemblies...")
    accession_list = selected_df["accession"].tolist()

    # Download in batches to avoid overlong command lines
    batch_size = 100
    all_extracted = []
    for i in range(0, len(accession_list), batch_size):
        batch = accession_list[i : i + batch_size]
        extracted = _download_accessions(batch, Path(output))
        all_extracted.extend(extracted)

    click.echo(f"Done! {len(all_extracted)} assemblies saved to {output}")


# -- accessions subcommand --

@main.command()
@click.argument("accessions_file", type=click.Path(exists=True, path_type=Path))
@click.option("--output", "-o", required=True, type=click.Path(path_type=Path),
              help="Output directory for downloaded assemblies.")
@click.option("--verbose", "-v", is_flag=True, default=False,
              help="Enable verbose logging.")
def accessions(accessions_file, output, verbose):
    """Fetch RefSeq assemblies for a list of accessions.

    ACCESSIONS_FILE is a text file with one accession per line.
    """
    _setup_logging(verbose)

    acc_list = [
        line.strip()
        for line in accessions_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    click.echo(f"Read {len(acc_list)} accessions from {accessions_file}")

    # Download in batches
    batch_size = 100
    all_extracted = []
    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i : i + batch_size]
        extracted = _download_accessions(batch, Path(output))
        all_extracted.extend(extracted)

    click.echo(f"Done! {len(all_extracted)} assemblies saved to {output}")


@main.command("list-accessions")
@click.argument("species_name")
@click.option("--verbose", "-v", is_flag=True, default=False,
              help="Enable verbose logging.")
def list_accessions(species_name, verbose):
    """List all RefSeq accessions for a species."""
    _setup_logging(verbose)

    metadata_df = _get_refseq_metadata(species_name)
    for acc in sorted(metadata_df["accession"].tolist()):
        click.echo(acc)


if __name__ == "__main__":
    main()
