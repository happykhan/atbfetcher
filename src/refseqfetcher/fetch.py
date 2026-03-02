"""Shared utilities for fetching from NCBI RefSeq.

This module wraps the NCBI ``datasets`` CLI tool to provide
programmatic access to RefSeq genome data.
"""

import json
import logging
import subprocess

logger = logging.getLogger(__name__)


def run_datasets(args: list[str], timeout: int = 300) -> str:
    """Run an NCBI datasets CLI command and return stdout.

    Parameters
    ----------
    args : list[str]
        Arguments to pass to the ``datasets`` command.
    timeout : int
        Command timeout in seconds.

    Returns
    -------
    str
        Standard output from the command.

    Raises
    ------
    RuntimeError
        If the command fails or ``datasets`` is not installed.
    """
    try:
        result = subprocess.run(
            ["datasets"] + args,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout,
        )
        return result.stdout
    except FileNotFoundError as exc:
        raise RuntimeError(
            "NCBI 'datasets' CLI not found. Install via pixi: pixi add ncbi-datasets-cli"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"datasets command failed: {exc.stderr}") from exc


def get_assembly_summary(taxon: str) -> list[dict]:
    """Get RefSeq assembly summaries for a taxon.

    Parameters
    ----------
    taxon : str
        Taxonomic name or NCBI taxon ID.

    Returns
    -------
    list[dict]
        List of assembly records from NCBI datasets.
    """
    output = run_datasets(
        [
            "summary",
            "genome",
            "taxon",
            taxon,
            "--assembly-source",
            "refseq",
            "--as-json-lines",
        ]
    )

    records = []
    for line in output.strip().splitlines():
        if not line.strip():
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            logger.warning("Skipping malformed JSON line")
            continue

    return records
