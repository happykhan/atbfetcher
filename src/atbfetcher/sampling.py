"""Stratified sampling of genomes by N50 and genome size.

Ensures the selected subset represents the full diversity of the
species by sampling proportionally across a 2D grid of Contig_N50
and Genome_Size quantile bins.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_N_BINS = 10


def stratified_sample(
    df: pd.DataFrame,
    n: int,
    seed: int = 42,
    n_bins: int = DEFAULT_N_BINS,
) -> pd.DataFrame:
    """Select n genomes using stratified sampling over N50 and genome size.

    Bins genomes into a quantile-based grid (default 10x10) of Contig_N50
    and Genome_Size, then samples proportionally from each bin.

    Parameters
    ----------
    df : pd.DataFrame
        Quality-filtered genomes. Must have ``Contig_N50`` and ``Genome_Size``.
    n : int
        Number of genomes to select.
    seed : int
        Random seed for reproducibility.
    n_bins : int
        Number of quantile bins per axis.

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

    # Create quantile-based bins for N50 and genome size
    df = df.copy()
    df["n50_bin"] = pd.qcut(
        df["Contig_N50"], q=n_bins, labels=False, duplicates="drop"
    )
    df["size_bin"] = pd.qcut(
        df["Genome_Size"], q=n_bins, labels=False, duplicates="drop"
    )

    # Group by the 2D grid
    grouped = df.groupby(["n50_bin", "size_bin"], observed=True)
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
        bin_df = grouped.get_group(bin_key)
        indices = rng.choice(len(bin_df), size=alloc, replace=False)
        selected_parts.append(bin_df.iloc[indices])

    result = pd.concat(selected_parts, ignore_index=True)

    # Drop temporary bin columns
    result = result.drop(columns=["n50_bin", "size_bin"], errors="ignore")

    logger.info(
        "Stratified sampling: selected %d from %d genomes across %d bins",
        len(result),
        len(df),
        len(allocations),
    )

    return result
