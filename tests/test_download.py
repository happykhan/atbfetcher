"""Tests for download and extraction logic."""

import tarfile
from pathlib import Path

import pandas as pd
import pytest

from atbfetcher.download import extract_samples, resolve_tarballs


class TestResolveTarballs:
    """Tests for mapping samples to tarballs."""

    def test_groups_by_tarball(self, sample_file_list_df):
        """Samples in the same tarball should be grouped together."""
        sample_ids = ["SAMN001", "SAMN002", "SAMN003"]
        result = resolve_tarballs(sample_ids, sample_file_list_df)

        assert "batch_001.tar.xz" in result
        assert "batch_002.tar.xz" in result

        # batch_001 should have SAMN001 and SAMN002
        batch1_samples = {e["sample"] for e in result["batch_001.tar.xz"]}
        assert batch1_samples == {"SAMN001", "SAMN002"}

        # batch_002 should have SAMN003
        batch2_samples = {e["sample"] for e in result["batch_002.tar.xz"]}
        assert batch2_samples == {"SAMN003"}

    def test_missing_samples_logged(self, sample_file_list_df):
        """Samples not in file list should be reported (not crash)."""
        sample_ids = ["SAMN001", "SAMN_MISSING"]
        result = resolve_tarballs(sample_ids, sample_file_list_df)

        # Should still find SAMN001
        all_samples = {
            e["sample"]
            for entries in result.values()
            for e in entries
        }
        assert "SAMN001" in all_samples
        assert "SAMN_MISSING" not in all_samples

    def test_empty_input_returns_empty(self, sample_file_list_df):
        """No sample IDs should return empty dict."""
        result = resolve_tarballs([], sample_file_list_df)
        assert result == {}


class TestExtractSamples:
    """Tests for extracting specific files from tar archives."""

    def test_extracts_requested_files(self, tmp_path):
        """Only requested files should be extracted."""
        # Create a test tar.xz with multiple files
        tar_path = tmp_path / "test.tar.xz"
        content_dir = tmp_path / "content"
        content_dir.mkdir()

        for name in ["file_a.fa.gz", "file_b.fa.gz", "file_c.fa.gz"]:
            (content_dir / name).write_text(f">{name}\nACGT\n")

        with tarfile.open(tar_path, "w:xz") as tar:
            for name in ["file_a.fa.gz", "file_b.fa.gz", "file_c.fa.gz"]:
                tar.add(content_dir / name, arcname=name)

        # Extract only two files
        output_dir = tmp_path / "output"
        result = extract_samples(tar_path, ["file_a.fa.gz", "file_c.fa.gz"], output_dir)

        assert len(result) == 2
        extracted_names = {p.name for p in result}
        assert extracted_names == {"file_a.fa.gz", "file_c.fa.gz"}

        # file_b should NOT be extracted
        assert not (output_dir / "file_b.fa.gz").exists()

    def test_missing_files_in_tarball(self, tmp_path):
        """Requesting a file not in the tarball should not crash."""
        tar_path = tmp_path / "test.tar.xz"
        content_dir = tmp_path / "content"
        content_dir.mkdir()

        (content_dir / "file_a.fa.gz").write_text(">seq\nACGT\n")
        with tarfile.open(tar_path, "w:xz") as tar:
            tar.add(content_dir / "file_a.fa.gz", arcname="file_a.fa.gz")

        output_dir = tmp_path / "output"
        result = extract_samples(
            tar_path, ["file_a.fa.gz", "file_missing.fa.gz"], output_dir
        )

        assert len(result) == 1
        assert result[0].name == "file_a.fa.gz"
