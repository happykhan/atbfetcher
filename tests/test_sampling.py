"""Tests for stratified sampling."""

import numpy as np
import pandas as pd
import pytest

from atbfetcher.sampling import stratified_sample


def _make_large_df(n: int = 500, seed: int = 0) -> pd.DataFrame:
    """Create a synthetic DataFrame for sampling tests."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "sample": [f"SAMN{i:06d}" for i in range(n)],
        "Contig_N50": rng.lognormal(mean=11, sigma=1.5, size=n).astype(int),
        "Genome_Size": rng.normal(5_000_000, 500_000, n).astype(int),
        "Completeness_Specific": rng.uniform(95, 100, n),
    })


class TestStratifiedSample:
    """Tests for the stratified sampling algorithm."""

    def test_returns_exact_count(self):
        """Should return exactly n samples when n < available."""
        df = _make_large_df(500)
        result = stratified_sample(df, n=100, seed=42)
        assert len(result) == 100

    def test_returns_all_when_fewer_than_n(self):
        """Should return all when fewer than n samples are available."""
        df = _make_large_df(50)
        result = stratified_sample(df, n=100, seed=42)
        assert len(result) == 50

    def test_reproducible_with_same_seed(self):
        """Same seed should produce identical results."""
        df = _make_large_df(500)
        result1 = stratified_sample(df, n=100, seed=42)
        result2 = stratified_sample(df, n=100, seed=42)
        pd.testing.assert_frame_equal(result1, result2)

    def test_different_seeds_differ(self):
        """Different seeds should produce different results."""
        df = _make_large_df(500)
        result1 = stratified_sample(df, n=100, seed=42)
        result2 = stratified_sample(df, n=100, seed=99)
        assert not result1["sample"].tolist() == result2["sample"].tolist()

    def test_covers_n50_range(self):
        """Selected samples should span the N50 range well."""
        df = _make_large_df(1000)
        result = stratified_sample(df, n=200, seed=42)

        # Use IQR-based check since N50 is lognormally distributed
        input_q25, input_q75 = df["Contig_N50"].quantile([0.25, 0.75])
        output_q25, output_q75 = result["Contig_N50"].quantile([0.25, 0.75])
        input_iqr = input_q75 - input_q25
        output_iqr = output_q75 - output_q25
        assert output_iqr >= 0.5 * input_iqr

    def test_covers_size_range(self):
        """Selected samples should span the genome size range."""
        df = _make_large_df(1000)
        result = stratified_sample(df, n=200, seed=42)

        input_range = df["Genome_Size"].max() - df["Genome_Size"].min()
        output_range = result["Genome_Size"].max() - result["Genome_Size"].min()
        assert output_range >= 0.7 * input_range

    def test_no_duplicate_samples(self):
        """Each sample should appear at most once."""
        df = _make_large_df(500)
        result = stratified_sample(df, n=100, seed=42)
        assert result["sample"].nunique() == len(result)

    def test_no_bin_columns_in_output(self):
        """Temporary bin columns should be dropped from output."""
        df = _make_large_df(500)
        result = stratified_sample(df, n=100, seed=42)
        assert "n50_bin" not in result.columns
        assert "size_bin" not in result.columns
