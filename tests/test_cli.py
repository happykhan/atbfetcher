"""Tests for the CLI interface."""

import pytest
from click.testing import CliRunner

from atbfetcher.cli import main


@pytest.fixture
def runner():
    """Click CLI test runner."""
    return CliRunner()


class TestCLI:
    """Integration tests for CLI commands."""

    def test_main_help(self, runner):
        """Main command should display help text."""
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "atbfetcher" in result.output

    def test_species_help(self, runner):
        """Species subcommand should display help text."""
        result = runner.invoke(main, ["species", "--help"])
        assert result.exit_code == 0
        assert "SPECIES_NAME" in result.output
        assert "--output" in result.output

    def test_mlst_help(self, runner):
        """MLST subcommand should display help text."""
        result = runner.invoke(main, ["mlst", "--help"])
        assert result.exit_code == 0
        assert "--scheme" in result.output

    def test_accessions_help(self, runner):
        """Accessions subcommand should display help text."""
        result = runner.invoke(main, ["accessions", "--help"])
        assert result.exit_code == 0
        assert "ACCESSIONS_FILE" in result.output

    def test_list_species_help(self, runner):
        """List-species subcommand should display help text."""
        result = runner.invoke(main, ["list-species", "--help"])
        assert result.exit_code == 0
        assert "--raw" in result.output

    def test_accessions_missing_file(self, runner):
        """Should error when accessions file doesn't exist."""
        result = runner.invoke(
            main, ["accessions", "/nonexistent/file.txt", "--output", "/tmp/out"]
        )
        assert result.exit_code != 0

    def test_version(self, runner):
        """Should display version info."""
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "version" in result.output.lower()
