# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-03-02

### Added
- `atbfetcher species` — stratified sampling by contig N50 and genome size
- `atbfetcher mlst` — MLST sequence type-based genome selection
- `atbfetcher accessions` — fetch specific genomes by accession ID
- `atbfetcher query` — filter genomes by species, country, year, host, isolation source via ATB SQLite metadata database
- `atbfetcher download-db` — download the ATB SQLite metadata database
- `atbfetcher list-species` — list available species with optional genome counts
- `atbfetcher species-count` — show HQ genome counts per species
- `atbfetcher list-countries` — list countries in the metadata database
- `atbfetcher list-hosts` — list host organisms in the metadata database
- `refseqfetcher species` — stratified sampling from NCBI RefSeq
- `refseqfetcher accessions` — fetch specific RefSeq assemblies
- Smart download source selection (AWS S3 vs OSF tarballs) with cost estimation
- Quality filtering via ATB HQ flag (default) or Qualibact per-species cutoffs
- Parquet-based metadata caching for fast reuse
- Multi-threaded XZ decompression for tarball downloads
- Rich CLI logging
