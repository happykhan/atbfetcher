# Download Source Selection

atbfetcher supports two download sources for AllTheBacteria assemblies:

| Source | Method | Best for |
|--------|--------|----------|
| **AWS** | Individual `.fa.gz` files from S3 | Small to medium batches |
| **OSF** | Extract from `tar.xz` archives | Large batches of a single species |

## How auto-selection works

When `--source auto` (the default), atbfetcher estimates the download time for
both sources and picks the faster one.

### Cost model

**AWS time** = `n_genomes × 0.4s`

Each genome is an individual HTTP request to the public S3 bucket
(`allthebacteria-assemblies.s3.eu-west-2.amazonaws.com`). Downloads run in
parallel (8 workers), so per-file time is amortised. Measured at ~0.35s/file
for 100 files.

**OSF time** = `n_tarballs × (60s + 25s)`

Each tarball (~3 GB) must be downloaded (~60s) and then decompressed + scanned
(~25s). Tarballs are processed sequentially because each is large.

### Breakeven point

Setting AWS = OSF:

```
n_genomes × 0.4 = n_tarballs × 85
n_genomes / n_tarballs = 212.5
```

OSF wins when there are **>213 genomes per tarball** on average. Because ATB
tarballs are grouped by species (closely related genomes are packed together),
this ratio improves rapidly as you fetch more genomes of the same species.

### Real-world examples

| Species | n | Tarballs | Genomes/tarball | AWS est. | OSF est. | Winner |
|---------|---|----------|-----------------|----------|----------|--------|
| S. aureus | 100 | 1 | 100 | 40s | 85s | AWS |
| S. aureus | 500 | 2 | 250 | 200s | 170s | OSF |
| S. aureus | 1000 | 2 | 500 | 400s | 170s | OSF |
| S. aureus | 5000 | 6 | 833 | 2000s | 510s | OSF |

### Measured benchmarks (100 S. aureus genomes)

| Source | Time | Notes |
|--------|------|-------|
| AWS | ~38s | 100 parallel downloads |
| OSF (cached) | ~10 min | 26 tarballs, extraction only |
| OSF (uncached) | ~47 min | 26 tarballs, download + extraction |

## CLI usage

```bash
# Auto-select (default) — estimates and picks the faster source
atbfetcher species "Staphylococcus aureus" -o output -n 100

# Force AWS (individual file downloads)
atbfetcher species "Staphylococcus aureus" -o output -n 100 --source aws

# Force OSF (tarball extraction)
atbfetcher species "Staphylococcus aureus" -o output -n 1000 --source osf
```

## Notes

- AWS downloads are `.fa.gz` (gzipped FASTA); tarball extractions are `.fa`
  (uncompressed FASTA from tar.xz)
- The auto-selection always assumes uncached tarballs for a fair comparison. If
  you have tarballs cached from a previous run and want to reuse them, use
  `--source osf`
- AWS S3 bucket is public and requires no authentication
