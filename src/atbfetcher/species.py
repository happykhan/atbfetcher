"""Species name handling for ATB genomes.

ATB uses GTDB taxonomy via Sylph, which appends suffixes like ``_A``, ``_B``
to species names for subspecies-level splits (e.g. "Enterobacter hormaechei_A").
This module provides utilities to clean those suffixes and match species
across different datasets (ATB, Qualibact, etc.).
"""

import re

import pandas as pd


def is_placeholder_species(name: str) -> bool:
    """Check if a species name is a GTDB unnamed placeholder.

    GTDB assigns numeric placeholder names like ``sp000746275`` to unnamed
    species. These are not useful for benchmarking and should be excluded.

    Parameters
    ----------
    name : str
        A species name to check.

    Returns
    -------
    bool
        True if the name looks like a GTDB placeholder (e.g. "Genus sp001234567").
    """
    return bool(re.search(r"\bsp\d{3,}", name))


def clean_species_name(name: str) -> str:
    """Strip GTDB subspecies suffixes from a species name.

    GTDB appends ``_A``, ``_B``, ``_C``, etc. to species names to denote
    subspecies-level genomic clusters. These suffixes are removed for display
    and cross-dataset matching.

    Parameters
    ----------
    name : str
        A species name, possibly with a GTDB suffix.

    Returns
    -------
    str
        The cleaned species name.

    Examples
    --------
    >>> clean_species_name("Enterobacter hormaechei_A")
    'Enterobacter hormaechei'
    >>> clean_species_name("Campylobacter_D jejuni")
    'Campylobacter jejuni'
    >>> clean_species_name("Escherichia coli")
    'Escherichia coli'
    """
    return re.sub(r"_[A-Z]\b", "", name)


def normalize_for_matching(name: str) -> str:
    """Normalize a species name for fuzzy cross-dataset matching.

    Strips GTDB suffixes and lowercases the name for comparison between
    ATB (GTDB taxonomy) and Qualibact (NCBI taxonomy).

    Parameters
    ----------
    name : str
        A species name from any source.

    Returns
    -------
    str
        Lowercased, suffix-stripped name.
    """
    return clean_species_name(name).lower().strip()


def list_species(species_calls_df: pd.DataFrame, raw: bool = False) -> list[str]:
    """Get a sorted list of unique species names from species calls.

    Parameters
    ----------
    species_calls_df : pd.DataFrame
        DataFrame with a ``species`` column (from Sylph species calls).
    raw : bool
        If True, return original GTDB names without cleaning.

    Returns
    -------
    list[str]
        Sorted unique species names.
    """
    col = "species"
    names = species_calls_df[col].dropna().unique()
    # Exclude GTDB placeholder species (e.g. "Genus sp000746275")
    names = [n for n in names if not is_placeholder_species(n)]

    if not raw:
        names = sorted(set(clean_species_name(n) for n in names))
    else:
        names = sorted(set(names))

    return names


def get_samples_for_species(
    name: str, species_calls_df: pd.DataFrame
) -> pd.DataFrame:
    """Get all sample IDs assigned to a species.

    Matches against both the raw GTDB name and the cleaned name,
    so "Enterobacter hormaechei" will match "Enterobacter hormaechei_A".

    Parameters
    ----------
    name : str
        Species name to search for.
    species_calls_df : pd.DataFrame
        DataFrame with ``sample`` and ``species`` columns.

    Returns
    -------
    pd.DataFrame
        Filtered rows matching the species.
    """
    if is_placeholder_species(name):
        return species_calls_df.iloc[0:0].copy()

    cleaned_input = normalize_for_matching(name)

    mask = species_calls_df["species"].apply(
        lambda x: (
            normalize_for_matching(str(x)) == cleaned_input
            and not is_placeholder_species(str(x))
        )
        if pd.notna(x)
        else False
    )
    return species_calls_df[mask].copy()
