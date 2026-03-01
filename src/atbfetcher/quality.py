"""Quality filtering using CheckM2 metrics and Qualibact cutoffs.

Applies per-species quality thresholds from Qualibact when available,
with sensible fallback defaults for species not in the Qualibact database.
"""

import logging

import pandas as pd

from atbfetcher.species import normalize_for_matching

logger = logging.getLogger(__name__)

# Default quality thresholds when species is not in Qualibact
DEFAULT_COMPLETENESS_MIN = 95.0
DEFAULT_CONTAMINATION_MAX = 5.0

# CheckM2 column names that map to Qualibact metric names.
# Both CheckM2 and Qualibact use the same metric names in the actual data.
CHECKM2_TO_QUALIBACT = {
    "Completeness_Specific": "Completeness_Specific",
    "Contamination": "Contamination",
    "Genome_Size": "Genome_Size",
    "GC_Content": "GC_Content",
}


def _find_qualibact_match(
    species_name: str, qualibact_cutoffs: dict
) -> dict | None:
    """Find matching Qualibact cutoffs for a species name.

    Normalizes both the query name and Qualibact species names for matching,
    stripping GTDB suffixes and lowercasing.

    Returns
    -------
    dict or None
        Metric cutoffs dict if found, None otherwise.
    """
    query = normalize_for_matching(species_name)
    for qb_species, cutoffs in qualibact_cutoffs.items():
        if normalize_for_matching(qb_species) == query:
            return cutoffs
    return None


def filter_by_quality(
    samples_df: pd.DataFrame,
    checkm2_df: pd.DataFrame,
    species_name: str,
    qualibact_cutoffs: dict,
) -> pd.DataFrame:
    """Filter samples by quality using CheckM2 metrics and Qualibact cutoffs.

    Merges the sample list with CheckM2 quality data, then applies either
    per-species Qualibact cutoffs or default thresholds.

    Parameters
    ----------
    samples_df : pd.DataFrame
        Samples to filter. Must have a ``sample`` column.
    checkm2_df : pd.DataFrame
        CheckM2 quality metrics. Must have ``sample`` plus quality columns.
    species_name : str
        Species name for looking up Qualibact cutoffs.
    qualibact_cutoffs : dict
        Nested dict from ``load_qualibact_cutoffs()``.

    Returns
    -------
    pd.DataFrame
        Quality-filtered samples with GC_Content and Genome_Size retained.
    """
    # Merge samples with CheckM2 quality data
    merged = samples_df.merge(checkm2_df, on="sample", how="inner")

    if merged.empty:
        logger.warning("No CheckM2 data found for any samples")
        return merged

    initial_count = len(merged)

    # Try to find species-specific cutoffs from Qualibact
    species_cutoffs = _find_qualibact_match(species_name, qualibact_cutoffs)

    if species_cutoffs is not None:
        logger.info("Using Qualibact cutoffs for %s", species_name)
        merged = _apply_qualibact_cutoffs(merged, species_cutoffs)
    else:
        logger.info(
            "No Qualibact cutoffs for %s, using defaults "
            "(completeness >= %.0f%%, contamination <= %.0f%%)",
            species_name,
            DEFAULT_COMPLETENESS_MIN,
            DEFAULT_CONTAMINATION_MAX,
        )
        merged = merged[
            (merged["Completeness_Specific"] >= DEFAULT_COMPLETENESS_MIN)
            & (merged["Contamination"] <= DEFAULT_CONTAMINATION_MAX)
        ]

    filtered_count = len(merged)
    logger.info(
        "Quality filter: %d -> %d samples (removed %d)",
        initial_count,
        filtered_count,
        initial_count - filtered_count,
    )

    return merged


def _apply_qualibact_cutoffs(
    df: pd.DataFrame, cutoffs: dict
) -> pd.DataFrame:
    """Apply Qualibact metric cutoffs to a DataFrame.

    For each metric in the cutoffs dict, filters rows where the corresponding
    CheckM2 column falls within [lower_bound, upper_bound].

    Contamination is handled specially: only the upper bound is applied.
    """
    result = df.copy()

    for checkm2_col, qualibact_metric in CHECKM2_TO_QUALIBACT.items():
        if qualibact_metric not in cutoffs:
            continue
        if checkm2_col not in result.columns:
            continue

        lower, upper = cutoffs[qualibact_metric]

        if qualibact_metric == "Contamination":
            # For contamination, only apply upper bound
            result = result[result[checkm2_col] <= upper]
        else:
            # Apply both bounds
            result = result[
                (result[checkm2_col] >= lower) & (result[checkm2_col] <= upper)
            ]

    return result
