"""Stratified sampling of genomes by assembly metrics.

Ensures the selected subset represents the full diversity of the
species by sampling proportionally across quantile bins. For ATB data,
bins on Contig_N50 and Genome_Size (2D grid). For RefSeq data (single
contig), bins on Genome_Size only (1D).
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_N_BINS = 10

# Default stratification axes for ATB (N50 + genome size) and RefSeq (genome size only)
ATB_COLUMNS = ("Contig_N50", "Genome_Size")
REFSEQ_COLUMNS = ("Genome_Size",)


def stratified_sample(
    df: pd.DataFrame,
    n: int,
    seed: int = 42,
    n_bins: int = DEFAULT_N_BINS,
    columns: tuple[str, ...] = ATB_COLUMNS,
) -> pd.DataFrame:
    """Select n genomes using stratified sampling over one or more metrics.

    Bins genomes into quantile bins along the given columns, then samples
    proportionally from each bin. Use ``ATB_COLUMNS`` (default) for N50 +
    genome size, or ``REFSEQ_COLUMNS`` for genome size only.

    Parameters
    ----------
    df : pd.DataFrame
        Quality-filtered genomes. Must have the columns listed in ``columns``.
    n : int
        Number of genomes to select.
    seed : int
        Random seed for reproducibility.
    n_bins : int
        Number of quantile bins per axis.
    columns : tuple[str, ...]
        Column names to stratify on. Default is ``("Contig_N50", "Genome_Size")``.
        For RefSeq, use ``("Genome_Size",)``.

    Returns
    -------
    pd.DataFrame
        Selected subset of genomes.
    """
    if len(df) <= n:
        logger.warning(
            "Requested %d samples but only %d available — returning all",
            n,
            len(df),
        )
        return df.copy()

    rng = np.random.default_rng(seed)

    # Create quantile-based bins for each stratification column
    df = df.copy()
    bin_col_names = []
    for col in columns:
        bin_name = f"_bin_{col}"
        bin_col_names.append(bin_name)
        df[bin_name] = pd.qcut(df[col], q=n_bins, labels=False, duplicates="drop")

    # Group by the bin grid (1D or 2D depending on columns)
    grouped = df.groupby(bin_col_names, observed=True)
    bin_counts = grouped.size()
    total = bin_counts.sum()

    # Calculate proportional allocation per bin (minimum 1 per non-empty bin)
    allocations = {}
    non_empty_bins = list(bin_counts.index)

    # First pass: proportional allocation
    for bin_key in non_empty_bins:
        proportion = bin_counts[bin_key] / total
        alloc = max(1, int(round(proportion * n)))
        # Don't allocate more than available in this bin
        alloc = min(alloc, bin_counts[bin_key])
        allocations[bin_key] = alloc

    # Adjust to hit exactly n
    current_total = sum(allocations.values())

    if current_total > n:
        # Reduce allocations from largest bins first
        sorted_bins = sorted(allocations, key=lambda k: allocations[k], reverse=True)
        for bin_key in sorted_bins:
            if current_total <= n:
                break
            reduction = min(allocations[bin_key] - 1, current_total - n)
            allocations[bin_key] -= reduction
            current_total -= reduction
    elif current_total < n:
        # Add to bins that have room, proportionally
        sorted_bins = sorted(
            allocations,
            key=lambda k: bin_counts[k] - allocations[k],
            reverse=True,
        )
        for bin_key in sorted_bins:
            if current_total >= n:
                break
            room = bin_counts[bin_key] - allocations[bin_key]
            addition = min(room, n - current_total)
            allocations[bin_key] += addition
            current_total += addition

    # Sample from each bin
    selected_parts = []
    for bin_key, alloc in allocations.items():
        # get_group needs a tuple for multi-column groupby, scalar for single
        group_key = (bin_key,) if not isinstance(bin_key, tuple) else bin_key
        bin_df = grouped.get_group(group_key)
        indices = rng.choice(len(bin_df), size=alloc, replace=False)
        selected_parts.append(bin_df.iloc[indices])

    result = pd.concat(selected_parts, ignore_index=True)

    # Drop temporary bin columns
    result = result.drop(columns=bin_col_names, errors="ignore")

    logger.info(
        "Stratified sampling: selected %d from %d genomes across %d bins",
        len(result),
        len(df),
        len(allocations),
    )

    return result
