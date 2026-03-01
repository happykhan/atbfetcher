"""Tests for quality filtering."""

import pandas as pd

from atbfetcher.quality import filter_by_quality


class TestFilterByQuality:
    """Tests for CheckM2 + Qualibact quality filtering."""

    def test_filter_with_qualibact_match(
        self, sample_species_calls_df, sample_checkm2_df, sample_qualibact_cutoffs
    ):
        """Species in Qualibact should use per-species cutoffs."""
        samples = sample_species_calls_df[
            sample_species_calls_df["species"] == "Escherichia coli"
        ]
        result = filter_by_quality(
            samples, sample_checkm2_df, "Escherichia coli", sample_qualibact_cutoffs
        )
        # SAMN001 (GC=50.5, size=5.1M, comp=99.5, cont=0.5) -> passes
        # SAMN002 (GC=50.8, size=4.9M, comp=97.2, cont=1.2) -> passes
        # SAMN007 (GC=50.2, size=5.05M, comp=96.0, cont=1.0) -> passes
        assert len(result) == 3
        assert set(result["sample"]) == {"SAMN001", "SAMN002", "SAMN007"}

    def test_filter_without_qualibact_match(
        self, sample_species_calls_df, sample_checkm2_df
    ):
        """Species not in Qualibact should use default thresholds."""
        samples = sample_species_calls_df[
            sample_species_calls_df["species"] == "Staphylococcus aureus"
        ]
        # No Qualibact entry → defaults: completeness >= 95, contamination <= 5
        result = filter_by_quality(
            samples, sample_checkm2_df, "Staphylococcus aureus", {}
        )
        # SAMN006 (comp=95.5, cont=2.1) -> passes defaults
        assert len(result) == 1
        assert result.iloc[0]["sample"] == "SAMN006"

    def test_filter_removes_low_quality(
        self, sample_species_calls_df, sample_checkm2_df
    ):
        """Low completeness / high contamination genomes should be filtered out."""
        samples = sample_species_calls_df[
            sample_species_calls_df["species"] == "Enterobacter hormaechei_A"
        ]
        # SAMN003 (comp=98.0, cont=0.8) -> passes defaults
        # SAMN004 (comp=60.0, cont=15.0) -> fails defaults
        result = filter_by_quality(
            samples, sample_checkm2_df, "Enterobacter hormaechei", {}
        )
        assert len(result) == 1
        assert result.iloc[0]["sample"] == "SAMN003"

    def test_empty_samples_returns_empty(self, sample_checkm2_df):
        """Empty input should return empty output."""
        empty = pd.DataFrame({"sample": []})
        result = filter_by_quality(empty, sample_checkm2_df, "Foo", {})
        assert len(result) == 0

    def test_qualibact_gc_filter(
        self, sample_species_calls_df, sample_checkm2_df, sample_qualibact_cutoffs
    ):
        """Qualibact GC_Content bounds should filter outliers."""
        # For E. coli, GC bounds are (49.0, 52.0)
        # All three E. coli samples have GC in range (50.2, 50.5, 50.8)
        samples = sample_species_calls_df[
            sample_species_calls_df["species"] == "Escherichia coli"
        ]
        result = filter_by_quality(
            samples, sample_checkm2_df, "Escherichia coli", sample_qualibact_cutoffs
        )
        for _, row in result.iterrows():
            assert 49.0 <= row["GC_Content"] <= 52.0
