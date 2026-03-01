"""Tests for the refseqfetcher CLI."""

from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from refseqfetcher.cli import main


@pytest.fixture
def runner():
    """Click CLI test runner."""
    return CliRunner()


class TestRefSeqFetcherCLI:
    """Integration tests for refseqfetcher CLI."""

    def test_main_help(self, runner):
        """Main command should display help text."""
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "refseqfetcher" in result.output

    def test_species_help(self, runner):
        """Species subcommand should display help text."""
        result = runner.invoke(main, ["species", "--help"])
        assert result.exit_code == 0
        assert "SPECIES_NAME" in result.output
        assert "--output" in result.output

    def test_accessions_help(self, runner):
        """Accessions subcommand should display help text."""
        result = runner.invoke(main, ["accessions", "--help"])
        assert result.exit_code == 0
        assert "ACCESSIONS_FILE" in result.output

    def test_accessions_missing_file(self, runner):
        """Should error when accessions file doesn't exist."""
        result = runner.invoke(main, [
            "accessions", "/nonexistent/file.txt", "--output", "/tmp/out"
        ])
        assert result.exit_code != 0

    def test_version(self, runner):
        """Should display version info."""
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0

    def test_list_accessions_help(self, runner):
        """List-accessions subcommand should display help text."""
        result = runner.invoke(main, ["list-accessions", "--help"])
        assert result.exit_code == 0
        assert "SPECIES_NAME" in result.output
