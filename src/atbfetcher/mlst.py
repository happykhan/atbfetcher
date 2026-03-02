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


def filter_by_mlst(
    quality_filtered_df: pd.DataFrame,
    mlst_df: pd.DataFrame,
    scheme: str | None,
    n: int,
    suspect_df: pd.DataFrame | None = None,
    seed: int = 42,
) -> pd.DataFrame:
    """Select genomes based on MLST sequence types.

    Ensures representation across STs. If n >= number of unique STs,
    selects one genome per ST (preferring highest completeness). If
    n < unique STs, randomly samples n STs and picks one genome each.

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

    Returns
    -------
    pd.DataFrame
        Selected genomes, one per ST.
    """
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

    logger.info("Unique resolved STs: %d, selecting up to %d", n_unique, n)

    selected_parts = []

    if n_unique > 0:
        if n >= n_unique:
            # One genome per ST — prefer highest completeness
            sts_to_use = unique_sts
        else:
            # Random subset of STs
            indices = rng.choice(n_unique, size=n, replace=False)
            sts_to_use = unique_sts[indices]

        for st in sts_to_use:
            st_samples = resolved[resolved["mlst_st"] == st]
            if "Completeness_Specific" in st_samples.columns:
                best = st_samples.sort_values("Completeness_Specific", ascending=False).iloc[0]
            else:
                best = st_samples.iloc[0]
            selected_parts.append(best)

    if not selected_parts:
        return scheme_df.head(0)

    result = pd.DataFrame(selected_parts)
    logger.info("MLST selection: %d genomes from %d STs", len(result), len(result))

    return result
