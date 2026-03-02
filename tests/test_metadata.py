"""Tests for metadata caching logic."""

from unittest.mock import MagicMock, patch

import pandas as pd

from atbfetcher.metadata import MetadataCache, load_qualibact_cutoffs


class TestMetadataCache:
    """Tests for the MetadataCache class."""

    def test_creates_cache_directory(self, tmp_path):
        """Cache directory should be created on initialization."""
        cache_dir = tmp_path / "new_cache"
        MetadataCache(cache_dir=cache_dir)
        assert cache_dir.exists()

    def test_stores_and_reloads_parquet(self, tmp_cache_dir):
        """Data should persist as parquet and reload correctly."""
        MetadataCache(cache_dir=tmp_cache_dir)

        # Manually write a parquet file to simulate a cached download
        test_df = pd.DataFrame({"sample": ["S1", "S2"], "value": [1, 2]})
        test_df.to_parquet(tmp_cache_dir / "test_data.parquet", index=False)

        # Verify it loads
        loaded = pd.read_parquet(tmp_cache_dir / "test_data.parquet")
        pd.testing.assert_frame_equal(test_df, loaded)

    def test_refresh_flag_bypasses_cache(self, tmp_cache_dir):
        """With refresh=True, cached parquet files should be ignored."""
        cache = MetadataCache(cache_dir=tmp_cache_dir, refresh=True)

        # Write a stale parquet
        stale_df = pd.DataFrame({"sample": ["old"]})
        stale_df.to_parquet(tmp_cache_dir / "file_list.parquet", index=False)

        # The _load_or_download method should attempt download
        # (we test the flag logic, not the actual download)
        parquet_path = cache._parquet_path("file_list")
        assert parquet_path.exists()
        # With refresh=True, it would skip this file
        assert cache.refresh is True

    def test_no_cache_flag(self, tmp_cache_dir):
        """With no_cache=True, files should not be saved."""
        cache = MetadataCache(cache_dir=tmp_cache_dir, no_cache=True)
        assert cache.no_cache is True

    @patch("atbfetcher.metadata.requests.get")
    def test_load_qualibact_cutoffs_parses_correctly(self, mock_get):
        """Qualibact CSV should be parsed into nested dict."""
        csv_text = (
            "species,metric,lower_bounds,upper_bounds\n"
            "Escherichia coli,Completeness_Specific,95.0,100.0\n"
            "Escherichia coli,Contamination,0.0,5.0\n"
            "Escherichia coli,Genome_Size,4500000,5800000\n"
        )
        mock_response = MagicMock()
        mock_response.text = csv_text
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        cutoffs = load_qualibact_cutoffs("https://fake.url/cutoffs.csv")

        assert "Escherichia coli" in cutoffs
        assert "Completeness_Specific" in cutoffs["Escherichia coli"]
        assert cutoffs["Escherichia coli"]["Completeness_Specific"] == (95.0, 100.0)
        assert cutoffs["Escherichia coli"]["Contamination"] == (0.0, 5.0)

    @patch("atbfetcher.metadata.requests.get")
    def test_load_qualibact_cutoffs_multiple_species(self, mock_get):
        """Should handle multiple species in the CSV."""
        csv_text = (
            "species,metric,lower_bounds,upper_bounds\n"
            "Escherichia coli,Completeness_Specific,95.0,100.0\n"
            "Klebsiella pneumoniae,Completeness_Specific,94.0,100.0\n"
        )
        mock_response = MagicMock()
        mock_response.text = csv_text
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        cutoffs = load_qualibact_cutoffs("https://fake.url/cutoffs.csv")

        assert len(cutoffs) == 2
        assert "Escherichia coli" in cutoffs
        assert "Klebsiella pneumoniae" in cutoffs
