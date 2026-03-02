"""Tests for download and extraction logic."""

import tarfile
from unittest.mock import MagicMock, patch

from atbfetcher.download import (
    AWS_BASE_URL,
    estimate_download_time,
    extract_samples,
    fetch_from_aws,
    resolve_tarballs,
)


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
        all_samples = {e["sample"] for entries in result.values() for e in entries}
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
        result = extract_samples(tar_path, ["file_a.fa.gz", "file_missing.fa.gz"], output_dir)

        assert len(result) == 1
        assert result[0].name == "file_a.fa.gz"


class TestFetchFromAws:
    """Tests for AWS S3 individual file downloads."""

    def test_downloads_files(self, tmp_path):
        """Should download files from AWS and return paths."""
        fake_content = b">seq\nACGT\n"

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.iter_content.return_value = [fake_content]
        mock_response.raise_for_status.return_value = None

        with patch("atbfetcher.download.requests.get", return_value=mock_response):
            result = fetch_from_aws(["SAMN001", "SAMN002"], tmp_path, max_workers=2)

        assert len(result) == 2
        names = {p.name for p in result}
        assert names == {"SAMN001.fa.gz", "SAMN002.fa.gz"}

    def test_skips_existing_files(self, tmp_path):
        """Should skip files that already exist."""
        existing = tmp_path / "SAMN001.fa.gz"
        existing.write_bytes(b">seq\nACGT\n")

        with patch("atbfetcher.download.requests.get") as mock_get:
            result = fetch_from_aws(["SAMN001"], tmp_path, max_workers=1)

        # Should not have made any HTTP requests
        mock_get.assert_not_called()
        assert len(result) == 1

    def test_handles_download_failure(self, tmp_path):
        """Failed downloads should return None and not crash."""
        import requests as req

        with patch(
            "atbfetcher.download.requests.get",
            side_effect=req.RequestException("Connection error"),
        ):
            result = fetch_from_aws(["SAMN_BAD"], tmp_path, max_workers=1)

        assert len(result) == 0
        assert not (tmp_path / "SAMN_BAD.fa.gz").exists()

    def test_correct_url_pattern(self, tmp_path):
        """Should use the correct AWS S3 URL pattern."""
        mock_response = MagicMock()
        mock_response.iter_content.return_value = [b"data"]
        mock_response.raise_for_status.return_value = None

        with patch("atbfetcher.download.requests.get", return_value=mock_response) as mock_get:
            fetch_from_aws(["SAMEA2445563"], tmp_path, max_workers=1)

        mock_get.assert_called_once_with(
            f"{AWS_BASE_URL}/SAMEA2445563.fa.gz",
            stream=True,
            timeout=120,
        )


class TestEstimateDownloadTime:
    """Tests for download time estimation logic."""

    def test_few_genomes_prefers_aws(self):
        """Small genome counts should prefer AWS."""
        method, aws_t, osf_t = estimate_download_time(n_genomes=10, n_tarballs=7)
        assert method == "aws"
        assert aws_t < osf_t

    def test_many_genomes_few_tarballs_prefers_osf(self):
        """When genomes are densely packed in few tarballs, OSF wins."""
        # 50000 genomes from 5 tarballs
        method, aws_t, osf_t = estimate_download_time(n_genomes=50000, n_tarballs=5)
        assert method == "osf"
        assert osf_t < aws_t

    def test_more_tarballs_increases_osf_time(self):
        """More tarballs should increase OSF estimated time."""
        _, _, few_t = estimate_download_time(n_genomes=50, n_tarballs=5)
        _, _, many_t = estimate_download_time(n_genomes=50, n_tarballs=20)
        assert many_t > few_t

    def test_returns_positive_times(self):
        """Estimated times should always be positive."""
        _, aws_t, osf_t = estimate_download_time(n_genomes=1, n_tarballs=1)
        assert aws_t > 0
        assert osf_t > 0
