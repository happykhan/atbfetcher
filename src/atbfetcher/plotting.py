"""Distribution plots for genome selection visualization.

Generates scatter plots of N50 vs genome size, highlighting the
selected subset against all high-quality genomes for a species.
"""

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for CLI use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

MAX_BACKGROUND_POINTS = 10_000


def plot_selection(
    all_hq_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    output_dir: Path,
    species_name: str,
) -> Path:
    """Create an N50 vs genome size scatter plot showing selected genomes.

    Grey dots show all high-quality genomes for the species. Colored dots
    highlight the selected subset. The plot is saved as a PNG file.

    Parameters
    ----------
    all_hq_df : pd.DataFrame
        All high-quality genomes (must have Contig_N50, Genome_Size).
    selected_df : pd.DataFrame
        Selected subset of genomes.
    output_dir : Path
        Directory to save the plot.
    species_name : str
        Species name for the plot title.

    Returns
    -------
    Path
        Path to the saved plot file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Subsample background points if too many
    bg_df = all_hq_df
    if len(bg_df) > MAX_BACKGROUND_POINTS:
        bg_df = bg_df.sample(n=MAX_BACKGROUND_POINTS, random_state=42)

    # Convert to kbp/Mbp for readability
    bg_size_mbp = bg_df["Genome_Size"] / 1_000_000
    sel_size_mbp = selected_df["Genome_Size"] / 1_000_000
    bg_n50_kbp = bg_df["Contig_N50"] / 1_000
    sel_n50_kbp = selected_df["Contig_N50"] / 1_000

    # Plot all HQ genomes as grey background
    ax.scatter(
        bg_n50_kbp,
        bg_size_mbp,
        c="lightgrey",
        s=8,
        alpha=0.5,
        label=f"All HQ (n={len(all_hq_df):,})",
        edgecolors="none",
    )

    # Plot selected genomes
    ax.scatter(
        sel_n50_kbp,
        sel_size_mbp,
        c="#2196F3",
        s=15,
        alpha=0.8,
        label=f"Selected (n={len(selected_df):,})",
        edgecolors="none",
    )

    ax.set_xlabel("Contig N50 (kbp)")
    ax.set_ylabel("Genome Size (Mbp)")
    ax.set_title(f"{species_name} — Genome Selection")
    ax.legend(loc="upper right")

    plt.tight_layout()

    safe_name = species_name.replace(" ", "_").replace("/", "_")
    plot_path = output_dir / f"{safe_name}_selection.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    logger.info("Saved selection plot to %s", plot_path)
    return plot_path
