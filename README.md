# atbfetcher

Fetch genomes from [AllTheBacteria](https://allthebacteria.org) (ATB) and [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) for constructing standardised benchmark datasets.

## Features

- **Species subsample**: Stratified sampling by contig N50 and genome size for representative subsets
- **MLST-based selection**: Pick genomes by MLST sequence type for phylogenetic diversity
- **Accession list**: Fetch specific genomes by accession ID
- **RefSeq mode**: Same workflow against NCBI RefSeq data (stratified by genome size)
- **Quality filtering**: Automatic quality control via [Qualibact](https://qualibact.org) per-species cutoffs and [CheckM2](https://github.com/chklovski/CheckM2) metrics
- **Caching**: Downloaded metadata cached locally as Parquet for fast reuse
- **Multi-threaded decompression**: Parallel XZ extraction for faster downloads
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

Requires `mlst_processed_all_samples.tsv.xz` in the `data/` directory.

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
| `--threads`, `-t` | CPU count - 1 | Threads for XZ decompression |

### `atbfetcher accessions`

Fetch specific assemblies by accession ID.

### `atbfetcher list-species`

Print available species names. Use `--raw` for original GTDB names, `--count` to show genome counts ordered highest to lowest.

### `atbfetcher species-count`

Show the number of HQ genomes per species. Use `--top N` to limit output.

### `refseqfetcher species`

Fetch a stratified subsample from NCBI RefSeq (stratified by genome size only since RefSeq genomes are typically single-contig).

### `refseqfetcher accessions`

Fetch specific RefSeq assemblies by accession.

## How It Works

### Data Sources

| File | Source | Purpose |
|------|--------|---------|
| `file_list.all.latest.tsv.gz` | OSF | Maps samples to download tarballs |
| `species_calls.tsv.gz` | OSF | Sylph species assignments (GTDB taxonomy) |
| `checkm2.tsv.gz` | OSF | Quality metrics (completeness, contamination, etc.) |
| Qualibact cutoffs | qualibact.org | Per-species quality thresholds |
| `mlst_processed_all_samples.tsv.xz` | User-provided | MLST assignments for all samples |

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
    download.py       # Tarball download/extraction (multi-threaded XZ)
    mlst.py           # MLST-based selection
  refseqfetcher/      # RefSeq genome fetching
    cli.py            # RefSeq CLI
    fetch.py          # NCBI datasets integration
tests/                # Pytest test suite
data/                 # Bundled reference data
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Write tests for your changes
4. Ensure all tests pass (`pixi run -e test test`)
5. Submit a pull request

## License

MIT License. See [LICENSE](LICENSE) for details.
