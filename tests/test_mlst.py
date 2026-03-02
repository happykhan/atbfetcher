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
        mlst_df = pd.DataFrame(
            {
                "sample": ["X1"],
                "mlst_scheme": ["-"],
                "mlst_st": ["-"],
            }
        )
        result = auto_detect_scheme(mlst_df, ["X1"])
        assert result is None

    def test_ignores_dash_scheme(self, sample_mlst_df):
        """Scheme '-' should be excluded from detection."""
        # Only pass samples that have no valid scheme
        empty_mlst = pd.DataFrame(
            {
                "sample": ["X1", "X2"],
                "mlst_scheme": ["-", "-"],
                "mlst_st": ["-", "-"],
            }
        )
        result = auto_detect_scheme(empty_mlst, ["X1", "X2"])
        assert result is None


class TestFilterByMlst:
    """Tests for MLST-based genome selection."""

    def _make_quality_df(self) -> pd.DataFrame:
        """Create a quality-filtered DataFrame for testing."""
        return pd.DataFrame(
            {
                "sample": ["SAMN001", "SAMN002", "SAMN007"],
                "species": ["Escherichia coli"] * 3,
                "Completeness_Specific": [99.5, 97.2, 96.0],
                "Contamination": [0.5, 1.2, 1.0],
                "GC_Content": [50.5, 50.8, 50.2],
                "Genome_Size": [5_100_000, 4_900_000, 5_050_000],
            }
        )

    def test_one_per_st(self, sample_mlst_df):
        """Should select one genome per resolved ST, excluding novel."""
        quality_df = self._make_quality_df()
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            seed=42,
        )
        # 2 resolved STs (10, 131); novel ST ("-") excluded
        assert len(result) == 2
        assert "-" not in result["mlst_st"].values

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
        suspect_df = pd.DataFrame(
            {
                "Suspect_MLST": ["ecoli_achtman_4(10)"],
                "Species_clash_with": ["Escherichia coli"],
            }
        )
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

    def test_excludes_novel_sts(self, sample_mlst_df):
        """Genomes with unresolved STs (st == '-') should be excluded."""
        quality_df = self._make_quality_df()
        result = filter_by_mlst(
            quality_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            seed=42,
        )
        assert "-" not in result["mlst_st"].values

    def test_empty_quality_df(self, sample_mlst_df):
        """Empty quality DataFrame should return empty result."""
        empty_df = pd.DataFrame(
            {
                "sample": pd.Series(dtype=str),
                "Completeness_Specific": pd.Series(dtype=float),
            }
        )
        result = filter_by_mlst(
            empty_df,
            sample_mlst_df,
            scheme="ecoli_achtman_4",
            n=10,
            seed=42,
        )
        assert len(result) == 0


# -- Larger fixture for strategy tests --


@pytest.fixture
def multi_st_data():
    """Fixture with multiple samples per ST to test strategy selection.

    STs and counts:
      ST131: 5 samples (most frequent)
      ST10:  3 samples
      ST73:  2 samples
      ST95:  1 sample
    Total: 11 samples, 4 unique STs
    """
    samples = []
    completeness = []
    sts = []

    # ST131: 5 samples
    for i in range(5):
        samples.append(f"S131_{i}")
        completeness.append(99.0 - i)
        sts.append("131")

    # ST10: 3 samples
    for i in range(3):
        samples.append(f"S10_{i}")
        completeness.append(98.0 - i)
        sts.append("10")

    # ST73: 2 samples
    for i in range(2):
        samples.append(f"S73_{i}")
        completeness.append(97.0 - i)
        sts.append("73")

    # ST95: 1 sample
    samples.append("S95_0")
    completeness.append(95.0)
    sts.append("95")

    quality_df = pd.DataFrame(
        {
            "sample": samples,
            "species": ["Escherichia coli"] * len(samples),
            "Completeness_Specific": completeness,
            "Contamination": [1.0] * len(samples),
        }
    )
    mlst_df = pd.DataFrame(
        {
            "sample": samples,
            "mlst_scheme": ["ecoli_achtman_4"] * len(samples),
            "mlst_st": sts,
            "mlst_status": ["EXACT"] * len(samples),
        }
    )
    return quality_df, mlst_df


class TestFrequencyStrategy:
    """Tests for the 'frequency' selection strategy."""

    def test_picks_most_common_first(self, multi_st_data):
        """Most frequent ST (131) should always appear in results."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=1, strategy="frequency"
        )
        assert len(result) == 1
        assert result.iloc[0]["mlst_st"] == "131"

    def test_one_per_st_when_n_equals_unique(self, multi_st_data):
        """With n=4 (unique STs), should pick one from each ST."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, strategy="frequency"
        )
        assert len(result) == 4
        assert result["mlst_st"].nunique() == 4

    def test_wraps_around(self, multi_st_data):
        """When n > unique STs, should pick multiple per ST in frequency order."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=6, strategy="frequency"
        )
        assert len(result) == 6
        # After one per ST (4), 2 more should go to the most frequent STs
        st_counts = result["mlst_st"].value_counts()
        assert st_counts["131"] == 2  # most frequent gets 2nd pick first

    def test_picks_by_completeness(self, multi_st_data):
        """Within each ST, should pick highest completeness first."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, strategy="frequency"
        )
        # The ST131 pick should be S131_0 (completeness 99.0)
        st131_rows = result[result["mlst_st"] == "131"]
        assert st131_rows.iloc[0]["sample"] == "S131_0"


class TestProportionalStrategy:
    """Tests for the 'proportional' selection strategy."""

    def test_common_sts_get_more_picks(self, multi_st_data):
        """ST131 (45% of genomes) should get more picks than ST95 (9%)."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=8, strategy="proportional"
        )
        assert len(result) == 8
        st_counts = result["mlst_st"].value_counts()
        assert st_counts["131"] > st_counts["95"]

    def test_minimum_one_per_st(self, multi_st_data):
        """When n >= unique STs, every ST should get at least 1 pick."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, strategy="proportional"
        )
        assert result["mlst_st"].nunique() == 4


class TestEqualStrategy:
    """Tests for the 'equal' selection strategy."""

    def test_equal_allocation(self, multi_st_data):
        """Each ST should get the same count (±1) when all have enough genomes."""
        quality_df, mlst_df = multi_st_data
        # n=4, 4 STs -> 1 each (all STs have >=1 sample)
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, strategy="equal"
        )
        assert len(result) == 4
        st_counts = result["mlst_st"].value_counts()
        assert st_counts.max() - st_counts.min() <= 1

    def test_remainder_to_top_frequency(self, multi_st_data):
        """With n=5 and 4 STs, remainder (1) goes to most frequent ST."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=5, strategy="equal"
        )
        assert len(result) == 5
        st_counts = result["mlst_st"].value_counts()
        # 5 // 4 = 1 each + 1 remainder to ST131
        assert st_counts["131"] == 2

    def test_redistribution_when_st_exhausted(self, multi_st_data):
        """If an ST can't fill its allocation, excess goes to others."""
        quality_df, mlst_df = multi_st_data
        # n=10: 10/4=2 each + 2 remainder -> ST131 gets 3, ST10 gets 3,
        # but ST95 only has 1 sample, so shortfall redistributed
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=10, strategy="equal"
        )
        assert len(result) == 10
        st_counts = result["mlst_st"].value_counts()
        assert st_counts["95"] == 1  # only 1 available


class TestRandomStrategy:
    """Tests for the 'random' selection strategy."""

    def test_reproducible_with_seed(self, multi_st_data):
        """Same seed should produce same result."""
        quality_df, mlst_df = multi_st_data
        r1 = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=3, seed=123, strategy="random"
        )
        r2 = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=3, seed=123, strategy="random"
        )
        assert list(r1["sample"]) == list(r2["sample"])

    def test_wraps_around(self, multi_st_data):
        """When n > unique STs, should pick extras from remaining genomes."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=6, strategy="random"
        )
        assert len(result) == 6
        # All 4 STs should be represented (one per ST + 2 extras)
        assert result["mlst_st"].nunique() == 4

    def test_n_less_than_unique_sts(self, multi_st_data):
        """When n < unique STs, should pick exactly n STs."""
        quality_df, mlst_df = multi_st_data
        result = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=2, seed=42, strategy="random"
        )
        assert len(result) == 2
        assert result["mlst_st"].nunique() == 2


class TestStrategyValidation:
    """Tests for strategy parameter validation."""

    def test_invalid_strategy_raises(self, multi_st_data):
        """Unknown strategy should raise ValueError."""
        quality_df, mlst_df = multi_st_data
        with pytest.raises(ValueError, match="Unknown strategy"):
            filter_by_mlst(
                quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, strategy="invalid"
            )

    def test_default_strategy_is_frequency(self, multi_st_data):
        """Default strategy should be frequency."""
        quality_df, mlst_df = multi_st_data
        # Call without explicit strategy (default)
        result_default = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, seed=42
        )
        result_explicit = filter_by_mlst(
            quality_df, mlst_df, scheme="ecoli_achtman_4", n=4, seed=42, strategy="frequency"
        )
        assert list(result_default["sample"]) == list(result_explicit["sample"])
