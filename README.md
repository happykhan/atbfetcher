# atbfetcher

[![Tests](https://github.com/happykhan/atbfetcher/actions/workflows/test.yml/badge.svg)](https://github.com/happykhan/atbfetcher/actions/workflows/test.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Fetch genomes from [AllTheBacteria](https://allthebacteria.org) (ATB) and [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) for constructing standardised benchmark datasets.

## Features

- **Species subsample**: Stratified sampling by contig N50 and genome size for representative subsets
- **MLST-based selection**: Pick genomes by MLST sequence type for phylogenetic diversity
- **Accession list**: Fetch specific genomes by accession ID
- **Metadata query**: Filter genomes by species, country, year, host, isolation source via ATB SQLite database
- **RefSeq mode**: Same workflow against NCBI RefSeq data (stratified by genome size)
- **Quality filtering**: Automatic quality control via ATB HQ flag (default) or [Qualibact](https://qualibact.org) per-species cutoffs
- **Smart download**: Auto-selects between AWS S3 (individual files) and OSF tarballs based on estimated speed
- **Caching**: Downloaded metadata cached locally as Parquet for fast reuse
- **Multi-threaded decompression**: Parallel XZ extraction for faster tarball downloads
- **Rich logging**: Colourful CLI output via [Rich](https://github.com/Textualize/rich)

## Installation

This project uses [Pixi](https://pixi.sh) to manage both conda packages (for `ncbi-datasets-cli`) and Python dependencies.

```bash
# Clone the repository
git clone https://github.com/happykhan/atbfetcher.git
cd atbfetcher

# Install all dependencies (creates isolated environment)
pixi install
```

### Docker

```bash
docker build -t atbfetcher .
docker run atbfetcher atbfetcher --help
```

## Quick Start

### List available species

```bash
pixi run atbfetcher list-species
pixi run atbfetcher list-species --raw     # show original GTDB names
pixi run atbfetcher list-species --count   # show genome counts, highest first
```

### Show species genome counts

```bash
pixi run atbfetcher species-count
pixi run atbfetcher species-count --top 20   # top 20 species
```

### Fetch a species subsample

```bash
pixi run atbfetcher species "Escherichia coli" \
  --output ./benchmark_ecoli \
  --n 1000 \
  --seed 42 \
  --threads 8
```

### Fetch by MLST sequence types

MLST data is auto-downloaded and cached on first use.

```bash
pixi run atbfetcher mlst "Escherichia coli" \
  --scheme ecoli_achtman_4 \
  --output ./benchmark_ecoli_mlst \
  --n 1000
```

### Fetch specific accessions

```bash
pixi run atbfetcher accessions accessions.txt \
  --output ./benchmark_custom
```

### Query by metadata (species, country, year, etc.)

The `query` command lets you filter genomes using the full ATB metadata database — including ENA metadata like country, collection date, host, and isolation source. This requires the ATB SQLite metadata database, which contains assembly info, quality metrics, Sylph species calls, and ENA sample metadata for all ATB genomes.

The database is ~2 GB to download (compressed) and ~27 GB uncompressed. Download it once with:

```bash
pixi run atbfetcher download-db
```

This downloads and decompresses the database to `~/.atbfetcher/`. You only need to do this once — the database is reused across all query commands.

```bash
# List E. coli from the UK collected 2020-2023
pixi run atbfetcher query --species "Escherichia coli" \
  --country "United Kingdom" --year-from 2020 --year-to 2023

# Fetch 500 S. aureus from blood samples
pixi run atbfetcher query --species "Staphylococcus aureus" \
  --isolation-source blood --n 500 --output ./saureus_blood

# Filter by CheckM2 quality metrics
pixi run atbfetcher query --species "Klebsiella pneumoniae" \
  --min-completeness 98 --max-contamination 2

# List available countries for a species
pixi run atbfetcher list-countries --species "Escherichia coli"

# List available host organisms
pixi run atbfetcher list-hosts --species "Staphylococcus aureus"
```

Without `--output`, the `query` command prints matching accessions to stdout (useful for piping or counting). With `--output`, it downloads the assemblies.

### Fetch from RefSeq

```bash
pixi run refseqfetcher species "Escherichia coli" \
  --output ./refseq_ecoli \
  --n 1000

pixi run refseqfetcher accessions accessions.txt \
  --output ./refseq_custom
```

## CLI Reference

### `atbfetcher species`

Fetch a stratified subsample of genomes for a species from ATB.

| Option | Default | Description |
|--------|---------|-------------|
| `--output`, `-o` | (required) | Output directory for assemblies |
| `--n`, `-n` | 1000 | Number of genomes to select |
| `--seed` | 42 | Random seed for reproducibility |
| `--source` | auto | Download source: `auto`, `aws`, or `osf` |
| `--threads`, `-t` | CPU count - 1 | Threads for XZ decompression |
| `--cache-dir` | `~/.atbfetcher` | Directory for cached metadata |
| `--no-cache` | False | Skip caching, download fresh |
| `--refresh` | False | Force re-download of cached metadata |
| `--verbose`, `-v` | False | Enable debug logging |

### `atbfetcher mlst`

Select genomes based on MLST sequence types, ensuring diversity across STs.

| Option | Default | Description |
|--------|---------|-------------|
| `--scheme` | auto-detect | MLST scheme name |
| `--output`, `-o` | (required) | Output directory |
| `--n`, `-n` | 1000 | Number of genomes to select |
| `--seed` | 42 | Random seed |
| `--strategy` | frequency | ST selection strategy: `frequency`, `proportional`, `equal`, or `random` |
| `--source` | auto | Download source: `auto`, `aws`, or `osf` |
| `--threads`, `-t` | CPU count - 1 | Threads for XZ decompression |

### `atbfetcher accessions`

Fetch specific assemblies by accession ID.

| Option | Default | Description |
|--------|---------|-------------|
| `--output`, `-o` | (required) | Output directory |
| `--source` | auto | Download source: `auto`, `aws`, or `osf` |
| `--threads`, `-t` | CPU count - 1 | Threads for XZ decompression |

### `atbfetcher query`

Query the ATB SQLite metadata database to select genomes by metadata filters. Without `--output`, prints matching accessions to stdout. With `--output`, downloads the assemblies.

| Option | Default | Description |
|--------|---------|-------------|
| `--species`, `-s` | None | Filter by species name |
| `--country`, `-c` | None | Filter by country (prefix match for `Country:Region` format) |
| `--year-from` | None | Minimum collection year (inclusive) |
| `--year-to` | None | Maximum collection year (inclusive) |
| `--host` | None | Filter by host organism (e.g. `Homo sapiens`) |
| `--isolation-source` | None | Filter by isolation source (e.g. `blood`, `stool`) |
| `--hq-only/--no-hq` | True | Only include ATB high-quality genomes |
| `--min-completeness` | None | Minimum CheckM2 completeness |
| `--max-contamination` | None | Maximum CheckM2 contamination |
| `--min-genome-size` | None | Minimum genome size in bases |
| `--max-genome-size` | None | Maximum genome size in bases |
| `--n`, `-n` | None | Maximum genomes to return (random subsample) |
| `--seed` | 42 | Random seed for reproducibility |
| `--output`, `-o` | None | Output directory (omit to just list accessions) |
| `--db-path` | `~/.atbfetcher/atb.metadata.202505.sqlite` | Path to ATB SQLite database |
| `--source` | auto | Download source: `auto`, `aws`, or `osf` |

### `atbfetcher list-countries`

List available countries in the ATB metadata database. Use `--species` to filter by species.

### `atbfetcher list-hosts`

List available host organisms in the ATB metadata database. Use `--species` to filter by species.

### `atbfetcher download-db`

Download the ATB SQLite metadata database (~2 GB download, ~27 GB uncompressed). Required for `query`, `list-countries`, and `list-hosts` commands. The database is stored in the cache directory and only needs to be downloaded once. Use `--refresh` to re-download.

### `atbfetcher list-species`

Print available species names. Use `--raw` for original GTDB names, `--count` to show genome counts ordered highest to lowest.

### `atbfetcher species-count`

Show the number of HQ genomes per species. Use `--top N` to limit output.

### `refseqfetcher species`

Fetch a stratified subsample from NCBI RefSeq (stratified by genome size only since RefSeq genomes are typically single-contig).

### `refseqfetcher accessions`

Fetch specific RefSeq assemblies by accession.

## How It Works

### Download Sources

atbfetcher can download assemblies from two sources:

- **AWS S3** — individual `.fa.gz` files from a public bucket (fast for small batches)
- **OSF tarballs** — `tar.xz` archives containing many assemblies (efficient for large batches)

By default (`--source auto`), atbfetcher estimates the download time for both and picks the faster option. ATB tarballs are grouped by species, so for large species-specific fetches (>500 genomes), tarballs are often faster. See [docs/download_sources.md](docs/download_sources.md) for details on the cost model.

### Data Sources

| File | Source | Purpose |
|------|--------|---------|
| `file_list.all.latest.tsv.gz` | OSF | Maps samples to download tarballs |
| `species_calls.tsv.gz` | OSF | Sylph species assignments (GTDB taxonomy) |
| `checkm2.tsv.gz` | OSF | Quality metrics (completeness, contamination, etc.) |
| Qualibact cutoffs | qualibact.org | Per-species quality thresholds |
| `mlst_processed_all_samples.tsv.xz` | Cloudflare R2 | MLST assignments for all samples (auto-downloaded) |
| `atb.metadata.202505.sqlite.xz` | [OSF](https://osf.io/download/my56u/) | SQLite database with assembly, ENA, CheckM2, and Sylph metadata |

### ATB Metadata Database

The `query` command uses the [ATB SQLite metadata database](https://allthebacteria.readthedocs.io/en/latest/metadata_sqlite.html), a comprehensive ~27 GB database containing assembly metadata, CheckM2 quality metrics, Sylph species calls, and ENA sample metadata (country, collection date, host, isolation source, etc.) for all ATB genomes. Run `atbfetcher download-db` to download it.

The other commands (`species`, `mlst`, `accessions`) use lighter-weight TSV metadata files (~50-200 MB each) that are auto-downloaded and cached as Parquet on first use — no manual setup required.

### Species Naming

ATB uses [GTDB taxonomy](https://gtdb.ecogenomic.org/) via [Sylph](https://github.com/bluenote-1577/sylph). GTDB appends suffixes like `_A`, `_B` to species names for subspecies-level genomic clusters. These suffixes are automatically stripped for display and matching. GTDB placeholder species (e.g. `sp000746275`) are excluded from all outputs.

### Quality Filtering

1. If the species has [Qualibact](https://qualibact.org) cutoffs, those per-species thresholds are applied
2. Otherwise, defaults are used: completeness >= 95%, contamination <= 5%

### Stratified Sampling

**ATB mode**: Genomes are binned into a quantile-based 10x10 grid of contig N50 and genome size. Samples are drawn proportionally from each bin to ensure the selected subset represents the full assembly quality and size diversity.

**RefSeq mode**: Since RefSeq genomes are typically single-contig, sampling uses genome size only (1D binning) rather than N50.

## Development

### Running Tests

```bash
pixi run -e test test
```

### Project Structure

```
src/
  atbfetcher/         # ATB genome fetching
    cli.py            # Click CLI with subcommands
    metadata.py       # Download/cache metadata files
    species.py        # Species name cleaning & filtering
    quality.py        # Quality filtering
    sampling.py       # Stratified sampling (N50+size or size-only)
    plotting.py       # N50 vs genome size plots
    download.py       # AWS S3 + tarball download/extraction
    mlst.py           # MLST-based selection
    query.py          # SQLite metadata query
  refseqfetcher/      # RefSeq genome fetching
    cli.py            # RefSeq CLI
    fetch.py          # NCBI datasets integration
tests/                # Pytest test suite
data/                 # Bundled reference data
docs/                 # Documentation
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Write tests for your changes
4. Ensure all tests pass (`pixi run -e test test`)
5. Submit a pull request

## Citation

If you use AllTheBacteria data, please cite the preprint:

> Hunt M et al. AllTheBacteria – all bacterial genomes assembled, available, and searchable. *bioRxiv* (2025). [https://doi.org/10.1101/2024.03.08.584059](https://doi.org/10.1101/2024.03.08.584059)

## Acknowledgements

- **Martin Hunt** & **Zamin Iqbal** for creating and maintaining AllTheBacteria
- **Jake Lacey** & **Torsten Seemann** for providing MLST typing information used in this tool

## License

MIT License. See [LICENSE](LICENSE) for details.
