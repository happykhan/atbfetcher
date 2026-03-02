"""MLST-based genome selection.

Selects genomes based on Multi-Locus Sequence Typing (MLST), ensuring
representation across sequence types (STs). Handles suspect scheme
contaminations and supports auto-detection of the most common MLST
scheme for a species.
"""

import logging
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

STRATEGIES = ("frequency", "proportional", "equal", "random")


def load_suspect_contaminations(path: Path | None = None) -> pd.DataFrame:
    """Load the suspect MLST scheme contaminations file.

    Parameters
    ----------
    path : Path, optional
        Path to the CSV file. If None, loads the bundled data file.

    Returns
    -------
    pd.DataFrame
        Columns: Suspect_MLST, Species_clash_with.
    """
    if path is None:
        data_dir = resources.files("atbfetcher").parent.parent / "data"
        path = data_dir / "Suspect_scheme_contaminations.csv"

    return pd.read_csv(path)


def auto_detect_scheme(mlst_df: pd.DataFrame, sample_ids: list[str]) -> str | None:
    """Detect the most common MLST scheme for a set of samples.

    Parameters
    ----------
    mlst_df : pd.DataFrame
        MLST results with columns: sample, mlst_scheme.
    sample_ids : list[str]
        Sample IDs to consider.

    Returns
    -------
    str or None
        Most common scheme name, or None if no schemes found.
    """
    matched = mlst_df[(mlst_df["sample"].isin(sample_ids)) & (mlst_df["mlst_scheme"] != "-")]
    if matched.empty:
        return None

    scheme_counts = matched["mlst_scheme"].value_counts()
    best_scheme = scheme_counts.index[0]
    logger.info(
        "Auto-detected MLST scheme: %s (%d samples)",
        best_scheme,
        scheme_counts.iloc[0],
    )
    return best_scheme


def _remove_suspect_combinations(
    df: pd.DataFrame,
    species_name: str,
    suspect_df: pd.DataFrame,
) -> pd.DataFrame:
    """Remove samples whose MLST scheme+ST is flagged as a suspect contamination.

    Parameters
    ----------
    df : pd.DataFrame
        MLST-annotated samples with mlst_scheme and mlst_st columns.
    species_name : str
        Current species name.
    suspect_df : pd.DataFrame
        Suspect contaminations (Suspect_MLST, Species_clash_with).

    Returns
    -------
    pd.DataFrame
        Samples with suspect entries removed.
    """
    if suspect_df.empty:
        return df

    # Build set of suspect MLST scheme(ST) patterns for this species
    suspect_patterns = set()
    for _, row in suspect_df.iterrows():
        clash_species = str(row["Species_clash_with"]).lower()
        if clash_species in species_name.lower():
            suspect_patterns.add(row["Suspect_MLST"])

    if not suspect_patterns:
        return df

    # Parse scheme(ST) from suspect patterns, e.g. "aeromonas(2363)"
    suspect_scheme_st = set()
    for pattern in suspect_patterns:
        if "(" in pattern and ")" in pattern:
            scheme = pattern.split("(")[0]
            st = pattern.split("(")[1].rstrip(")")
            suspect_scheme_st.add((scheme, st))

    initial = len(df)
    mask = df.apply(
        lambda row: (str(row["mlst_scheme"]), str(row["mlst_st"])) not in suspect_scheme_st,
        axis=1,
    )
    result = df[mask]

    removed = initial - len(result)
    if removed > 0:
        logger.info("Removed %d suspect MLST combinations", removed)

    return result


def _sort_by_completeness(df: pd.DataFrame) -> pd.DataFrame:
    """Sort a DataFrame by descending completeness if the column exists."""
    if "Completeness_Specific" in df.columns:
        return df.sort_values("Completeness_Specific", ascending=False)
    return df


def _select_by_frequency(
    resolved: pd.DataFrame, n: int, rng: np.random.Generator
) -> list[pd.Series]:
    """Round-robin by ST frequency (most common first), wrapping as needed.

    Each pass iterates through STs in descending frequency order and picks
    the next-best genome (by completeness) that hasn't been picked yet.
    """
    st_counts = resolved["mlst_st"].value_counts()  # descending by default
    st_order = st_counts.index.tolist()

    selected: list[pd.Series] = []
    picked_indices: set = set()

    while len(selected) < n:
        added_this_round = False
        for st in st_order:
            if len(selected) >= n:
                break
            st_samples = _sort_by_completeness(resolved[resolved["mlst_st"] == st])
            for idx, row in st_samples.iterrows():
                if idx not in picked_indices:
                    selected.append(row)
                    picked_indices.add(idx)
                    added_this_round = True
                    break
        if not added_this_round:
            break  # exhausted all genomes

    return selected


def _select_proportional(
    resolved: pd.DataFrame, n: int, rng: np.random.Generator
) -> list[pd.Series]:
    """Allocate picks per ST proportional to their frequency."""
    st_counts = resolved["mlst_st"].value_counts()
    total = st_counts.sum()
    n_unique = len(st_counts)

    # Compute raw proportional allocation
    allocations: dict[str, int] = {}
    for st, count in st_counts.items():
        allocations[st] = (
            max(1, round(count / total * n)) if n >= n_unique else round(count / total * n)
        )

    # Adjust to hit exactly n
    current_total = sum(allocations.values())
    # Sort STs by frequency descending for consistent tie-breaking
    sts_by_freq = st_counts.index.tolist()

    if current_total > n:
        # Remove from smallest STs first (reverse order)
        for st in reversed(sts_by_freq):
            if current_total <= n:
                break
            reduce = min(allocations[st] - (1 if n >= n_unique else 0), current_total - n)
            if reduce > 0:
                allocations[st] -= reduce
                current_total -= reduce
    elif current_total < n:
        # Add to largest STs first
        for st in sts_by_freq:
            if current_total >= n:
                break
            # Don't exceed available genomes for this ST
            available = len(resolved[resolved["mlst_st"] == st])
            can_add = min(available - allocations[st], n - current_total)
            if can_add > 0:
                allocations[st] += can_add
                current_total += can_add

    selected: list[pd.Series] = []
    for st in sts_by_freq:
        st_samples = _sort_by_completeness(resolved[resolved["mlst_st"] == st])
        count = min(allocations.get(st, 0), len(st_samples))
        for _, row in st_samples.head(count).iterrows():
            selected.append(row)

    return selected[:n]


def _select_equal(resolved: pd.DataFrame, n: int, rng: np.random.Generator) -> list[pd.Series]:
    """Equal genomes per ST, remainder distributed by frequency."""
    st_counts = resolved["mlst_st"].value_counts()
    sts_by_freq = st_counts.index.tolist()
    n_unique = len(sts_by_freq)

    if n_unique == 0:
        return []

    per_st = n // n_unique
    remainder = n % n_unique

    # STs that get an extra genome (top-frequency STs)
    extra_sts = set(sts_by_freq[:remainder])

    allocations: dict[str, int] = {}
    for st in sts_by_freq:
        allocations[st] = per_st + (1 if st in extra_sts else 0)

    # Redistribute if any ST has fewer genomes than allocated
    shortfall = 0
    for st in sts_by_freq:
        available = len(resolved[resolved["mlst_st"] == st])
        if available < allocations[st]:
            shortfall += allocations[st] - available
            allocations[st] = available

    # Distribute shortfall to STs that have room, in frequency order
    if shortfall > 0:
        for st in sts_by_freq:
            if shortfall <= 0:
                break
            available = len(resolved[resolved["mlst_st"] == st])
            room = available - allocations[st]
            if room > 0:
                give = min(room, shortfall)
                allocations[st] += give
                shortfall -= give

    selected: list[pd.Series] = []
    for st in sts_by_freq:
        st_samples = _sort_by_completeness(resolved[resolved["mlst_st"] == st])
        count = min(allocations.get(st, 0), len(st_samples))
        for _, row in st_samples.head(count).iterrows():
            selected.append(row)

    return selected[:n]


def _select_random(resolved: pd.DataFrame, n: int, rng: np.random.Generator) -> list[pd.Series]:
    """Uniform random ST selection; wraps around when n > unique STs."""
    unique_sts = resolved["mlst_st"].unique()
    n_unique = len(unique_sts)

    if n_unique == 0:
        return []

    selected: list[pd.Series] = []
    picked_indices: set = set()

    if n <= n_unique:
        # Random subset of STs, one genome per ST (highest completeness)
        indices = rng.choice(n_unique, size=n, replace=False)
        sts_to_use = unique_sts[indices]
        for st in sts_to_use:
            st_samples = _sort_by_completeness(resolved[resolved["mlst_st"] == st])
            row = st_samples.iloc[0]
            selected.append(row)
            picked_indices.add(st_samples.index[0])
    else:
        # One per ST (highest completeness), then random from remaining
        for st in unique_sts:
            st_samples = _sort_by_completeness(resolved[resolved["mlst_st"] == st])
            row = st_samples.iloc[0]
            selected.append(row)
            picked_indices.add(st_samples.index[0])

        # Remaining picks from unpicked genomes
        remaining_n = n - len(selected)
        if remaining_n > 0:
            remaining = resolved[~resolved.index.isin(picked_indices)]
            if len(remaining) > 0:
                pick_n = min(remaining_n, len(remaining))
                extra_indices = rng.choice(len(remaining), size=pick_n, replace=False)
                for idx in extra_indices:
                    selected.append(remaining.iloc[idx])

    return selected


def filter_by_mlst(
    quality_filtered_df: pd.DataFrame,
    mlst_df: pd.DataFrame,
    scheme: str | None,
    n: int,
    suspect_df: pd.DataFrame | None = None,
    seed: int = 42,
    strategy: str = "frequency",
) -> pd.DataFrame:
    """Select genomes based on MLST sequence types.

    Ensures representation across STs using a configurable selection strategy.
    Genomes with unresolved STs (novel or unknown, st == "-") are excluded.

    Parameters
    ----------
    quality_filtered_df : pd.DataFrame
        Quality-filtered samples with a ``sample`` column and quality metrics.
    mlst_df : pd.DataFrame
        MLST results (sample, mlst_scheme, mlst_st, mlst_status).
    scheme : str or None
        MLST scheme to filter by. If None, auto-detects most common.
    n : int
        Number of genomes to select.
    suspect_df : pd.DataFrame, optional
        Suspect contaminations to filter out.
    seed : int
        Random seed for reproducibility.
    strategy : str
        ST selection strategy: 'frequency', 'proportional', 'equal', or 'random'.

    Returns
    -------
    pd.DataFrame
        Selected genomes.
    """
    if strategy not in STRATEGIES:
        msg = f"Unknown strategy {strategy!r}, must be one of {STRATEGIES}"
        raise ValueError(msg)
    rng = np.random.default_rng(seed)

    # Join quality-filtered samples with MLST data
    merged = quality_filtered_df.merge(mlst_df, on="sample", how="inner")

    if merged.empty:
        logger.warning("No MLST data found for quality-filtered samples")
        return merged

    # Auto-detect scheme if not specified
    sample_ids = quality_filtered_df["sample"].tolist()
    if scheme is None:
        scheme = auto_detect_scheme(mlst_df, sample_ids)
        if scheme is None:
            logger.warning("Could not auto-detect MLST scheme")
            return merged.head(0)

    # Filter to specified scheme
    scheme_df = merged[merged["mlst_scheme"] == scheme].copy()
    logger.info("Samples with scheme %s: %d", scheme, len(scheme_df))

    if scheme_df.empty:
        logger.warning("No samples found for scheme %s", scheme)
        return scheme_df

    # Remove suspect combinations
    if suspect_df is not None:
        species_name = (
            quality_filtered_df.get("species", pd.Series()).iloc[0]
            if "species" in quality_filtered_df.columns
            else ""
        )
        scheme_df = _remove_suspect_combinations(scheme_df, species_name, suspect_df)

    # Exclude novel/unknown STs (st == "-")
    resolved = scheme_df[scheme_df["mlst_st"] != "-"]
    n_excluded = len(scheme_df) - len(resolved)
    if n_excluded > 0:
        logger.info("Excluded %d samples with unresolved STs", n_excluded)

    unique_sts = resolved["mlst_st"].unique()
    n_unique = len(unique_sts)

    logger.info("Unique resolved STs: %d, selecting up to %d (strategy: %s)", n_unique, n, strategy)

    dispatch = {
        "frequency": _select_by_frequency,
        "proportional": _select_proportional,
        "equal": _select_equal,
        "random": _select_random,
    }
    selected_parts = dispatch[strategy](resolved, n, rng)

    if not selected_parts:
        return scheme_df.head(0)

    result = pd.DataFrame(selected_parts)
    n_sts = result["mlst_st"].nunique()
    logger.info("MLST selection: %d genomes from %d STs", len(result), n_sts)

    return result
