"""Shared test fixtures for atbfetcher tests.

Provides small synthetic DataFrames mimicking real ATB metadata,
plus temporary directories and mock qualibact cutoffs.
"""

import pandas as pd
import pytest


@pytest.fixture
def sample_species_calls_df():
    """Small species calls DataFrame for testing.

    Includes an HQ column mimicking ATB's high-quality flag.
    SAMN004 and SAMN008 are not HQ.
    """
    return pd.DataFrame(
        {
            "sample": [
                "SAMN001",
                "SAMN002",
                "SAMN003",
                "SAMN004",
                "SAMN005",
                "SAMN006",
                "SAMN007",
                "SAMN008",
            ],
            "species": [
                "Escherichia coli",
                "Escherichia coli",
                "Enterobacter hormaechei_A",
                "Enterobacter hormaechei_A",
                "Klebsiella pneumoniae",
                "Staphylococcus aureus",
                "Escherichia coli",
                "Bacillus sp000746275",
            ],
            "HQ": ["T", "T", "T", "F", "T", "T", "T", "T"],
        }
    )


@pytest.fixture
def sample_checkm2_df():
    """Small CheckM2 quality DataFrame for testing."""
    return pd.DataFrame(
        {
            "sample": [
                "SAMN001",
                "SAMN002",
                "SAMN003",
                "SAMN004",
                "SAMN005",
                "SAMN006",
                "SAMN007",
            ],
            "Completeness_Specific": [99.5, 97.2, 98.0, 60.0, 99.1, 95.5, 96.0],
            "Contamination": [0.5, 1.2, 0.8, 15.0, 0.3, 2.1, 1.0],
            "Genome_Size": [
                5_100_000,
                4_900_000,
                4_800_000,
                5_200_000,
                5_600_000,
                2_800_000,
                5_050_000,
            ],
            "GC_Content": [50.5, 50.8, 55.1, 54.9, 57.3, 32.8, 50.2],
            "Contig_N50": [
                150_000,
                120_000,
                95_000,
                10_000,
                200_000,
                80_000,
                130_000,
            ],
        }
    )


@pytest.fixture
def sample_file_list_df():
    """Small file list DataFrame mapping samples to tarballs."""
    return pd.DataFrame(
        {
            "sample": ["SAMN001", "SAMN002", "SAMN003", "SAMN005", "SAMN007"],
            "filename_in_tar_xz": [
                "SAMN001.fa.gz",
                "SAMN002.fa.gz",
                "SAMN003.fa.gz",
                "SAMN005.fa.gz",
                "SAMN007.fa.gz",
            ],
            "tar_xz": [
                "batch_001.tar.xz",
                "batch_001.tar.xz",
                "batch_002.tar.xz",
                "batch_002.tar.xz",
                "batch_001.tar.xz",
            ],
            "tar_xz_url": [
                "https://example.com/batch_001.tar.xz",
                "https://example.com/batch_001.tar.xz",
                "https://example.com/batch_002.tar.xz",
                "https://example.com/batch_002.tar.xz",
                "https://example.com/batch_001.tar.xz",
            ],
            "tar_xz_md5": [
                "abc123",
                "abc123",
                "def456",
                "def456",
                "abc123",
            ],
        }
    )


@pytest.fixture
def sample_qualibact_cutoffs():
    """Qualibact cutoffs dict for a few species."""
    return {
        "Escherichia coli": {
            "Completeness_Specific": (95.0, 100.0),
            "Contamination": (0.0, 5.0),
            "Genome_Size": (4_500_000, 5_800_000),
            "GC_Content": (49.0, 52.0),
        },
        "Klebsiella pneumoniae": {
            "Completeness_Specific": (95.0, 100.0),
            "Contamination": (0.0, 5.0),
            "Genome_Size": (5_000_000, 6_500_000),
            "GC_Content": (56.0, 58.5),
        },
    }


@pytest.fixture
def sample_mlst_df():
    """Small MLST DataFrame for testing."""
    return pd.DataFrame(
        {
            "sample": [
                "SAMN001",
                "SAMN002",
                "SAMN003",
                "SAMN005",
                "SAMN006",
                "SAMN007",
            ],
            "mlst_scheme": [
                "ecoli_achtman_4",
                "ecoli_achtman_4",
                "ehormaechei",
                "klebsiella",
                "saureus",
                "ecoli_achtman_4",
            ],
            "mlst_st": ["10", "131", "15", "258", "8", "-"],
            "mlst_status": [
                "EXACT",
                "EXACT",
                "EXACT",
                "EXACT",
                "EXACT",
                "NOVEL",
            ],
        }
    )


@pytest.fixture
def tmp_cache_dir(tmp_path):
    """Temporary directory for cache tests."""
    cache = tmp_path / "cache"
    cache.mkdir()
    return cache
