"""Terminal-aware tabular output formatting.

Detects whether stdout is a TTY and auto-selects between pretty table
output (for humans) and TSV (for pipes/scripts).

Supported formats: ``table``, ``tsv``, ``csv``, ``json``, ``auto``.

``auto`` resolves to ``table`` when stdout is a TTY, otherwise ``tsv``.
"""

from __future__ import annotations

import sys

import pandas as pd

FORMATS = ("table", "tsv", "csv", "json", "auto")


def _is_tty() -> bool:
    """Return True when stdout is connected to a terminal."""
    return sys.stdout.isatty()


def resolve_format(fmt: str) -> str:
    """Resolve 'auto' to a concrete format based on terminal detection.

    Parameters
    ----------
    fmt : str
        One of ``table|tsv|csv|json|auto``.

    Returns
    -------
    str
        A concrete format string (never ``auto``).
    """
    if fmt == "auto":
        return "table" if _is_tty() else "tsv"
    return fmt


def print_dataframe(df: pd.DataFrame, fmt: str = "auto") -> None:
    """Print a DataFrame to stdout in the requested format.

    Parameters
    ----------
    df : pd.DataFrame
        Data to print.
    fmt : str
        Output format: ``table|tsv|csv|json|auto``.
    """
    resolved = resolve_format(fmt)

    if resolved == "table":
        _print_table(df)
    elif resolved == "tsv":
        print(df.to_csv(sep="\t", index=False), end="")
    elif resolved == "csv":
        print(df.to_csv(index=False), end="")
    elif resolved == "json":
        print(df.to_json(orient="records", indent=2))
    else:
        # Unknown format — fall back to TSV
        print(df.to_csv(sep="\t", index=False), end="")


def _print_table(df: pd.DataFrame) -> None:
    """Print a DataFrame as a plain-text aligned table.

    Uses pandas built-in to_string for simplicity (no extra deps).
    """
    if df.empty:
        print("(no results)")
        return
    print(df.to_string(index=False))
