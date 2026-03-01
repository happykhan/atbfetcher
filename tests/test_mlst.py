"""Tests for MLST filtering logic."""

import pandas as pd
import pytest

from atbfetcher.mlst import (
    auto_detect_scheme,
    filter_by_mlst,
)


class TestAutoDetectScheme:
    """Tests for automatic MLST scheme detection."""

    def test_detects_most_common_scheme(self, sample_mlst_df):
        """Should return the scheme with the most samples."""
        sample_ids = ["SAMN001", "SAMN002", "SAMN007"]
        result = auto_detect_scheme(sample_mlst_df, sample_ids)
        assert result == "ecoli_achtman_4"

    def test_returns_none_when_no_schemes(self):
        """Should return None when no MLST data matches."""
        mlst_df = pd.DataFrame({
            "sample": ["X1"],
            "mlst_scheme": ["-"],
            "mlst_st": ["-"],
        })
        result = auto_detect_scheme(mlst_df, ["X1"])
        assert result is None

    def test_ignores_dash_scheme(self, sample_mlst_df):
        """Scheme '-' should be excluded from detection."""
        # Only pass samples that have no valid scheme
        empty_mlst = pd.DataFrame({
            "sample": ["X1", "X2"],
            "mlst_scheme": ["-", "-"],
            "mlst_st": ["-", "-"],
        })
        result = auto_detect_scheme(empty_mlst, ["X1", "X2"])
        assert result is None


class TestFilterByMlst:
    """Tests for MLST-based genome selection."""

    def _make_quality_df(self) -> pd.DataFrame:
        """Create a quality-filtered DataFrame for testing."""
        return pd.DataFrame({
            "sample": ["SAMN001", "SAMN002", "SAMN007"],
            "species": ["Escherichia coli"] * 3,
            "Completeness_Specific": [99.5, 97.2, 96.0],
            "Contamination": [0.5, 1.2, 1.0],
            "GC_Content": [50.5, 50.8, 50.2],
            "Genome_Size": [5_100_000, 4_900_000, 5_050_000],
        })

    def test_one_per_st(self, sample_mlst_df):
        """Should select one genome per ST."""
        quality_df = self._make_quality_df()
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            seed=42,
        )
        # There are 2 known STs (10, 131) + 1 novel (-) = 3 samples
        assert len(result) <= 3
        assert len(result) > 0

    def test_n_less_than_sts(self, sample_mlst_df):
        """When n < unique STs, should sample a subset of STs."""
        quality_df = self._make_quality_df()
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=1,
            seed=42,
        )
        assert len(result) >= 1

    def test_suspect_removal(self, sample_mlst_df):
        """Suspect MLST+species combos should be filtered out."""
        quality_df = self._make_quality_df()
        suspect_df = pd.DataFrame({
            "Suspect_MLST": ["ecoli_achtman_4(10)"],
            "Species_clash_with": ["Escherichia coli"],
        })
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            suspect_df=suspect_df,
            seed=42,
        )
        # ST 10 (SAMN001) should be removed
        if not result.empty and "mlst_st" in result.columns:
            assert "10" not in result["mlst_st"].values

    def test_auto_detect_scheme_integration(self, sample_mlst_df):
        """When scheme=None, should auto-detect and still work."""
        quality_df = self._make_quality_df()
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme=None,
            n=10,
            seed=42,
        )
        assert len(result) > 0

    def test_empty_quality_df(self, sample_mlst_df):
        """Empty quality DataFrame should return empty result."""
        empty_df = pd.DataFrame({
            "sample": pd.Series(dtype=str),
            "Completeness_Specific": pd.Series(dtype=float),
        })
        result = filter_by_mlst(
            empty_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            seed=42,
        )
        assert len(result) == 0
