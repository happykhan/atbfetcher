"""Tests for new CLI features:
1. Species typo suggestions
2. Terminal-aware output formatting
3. info subcommand
4. Dry-run for downloads
5. mlst-query subcommand
6. TOML filter files
7. summarise subcommand
8. Config file management
"""

import sqlite3

import pandas as pd
import pytest
from click.testing import CliRunner

from atbfetcher.cli import main, _suggest_species, _load_toml_filters
from atbfetcher.config import (
    CONFIG_PATH,
    default_config,
    load_config,
    write_config,
)
from atbfetcher.output import print_dataframe, resolve_format


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def test_db(tmp_path):
    """Reusable small test SQLite database."""
    db_path = tmp_path / "test.sqlite"
    conn = sqlite3.connect(str(db_path))

    conn.execute("""
        CREATE TABLE assembly (
            sample_accession TEXT PRIMARY KEY,
            run_accession TEXT,
            assembly_accession TEXT,
            asm_fasta_on_osf INTEGER,
            sylph_species TEXT,
            hq_filter TEXT,
            osf_tarball_filename TEXT,
            osf_tarball_url TEXT,
            aws_url TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE run (
            run_accession TEXT PRIMARY KEY,
            sample_accession TEXT,
            pass INTEGER
        )
    """)
    conn.execute("""
        CREATE TABLE ena_202505_used (
            run_accession TEXT,
            sample_accession TEXT,
            country TEXT,
            collection_date TEXT,
            host TEXT,
            isolation_source TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE checkm2 (
            sample_accession TEXT PRIMARY KEY,
            Completeness_Specific REAL,
            Contamination REAL,
            Genome_Size INT,
            GC_Content FLOAT,
            Contig_N50 INT
        )
    """)

    assemblies = [
        ("SAMN001", "SRR001", "GCA_001", 1, "Escherichia coli", "PASS",
         "batch_001.tar.xz", "https://osf.io/batch_001", "https://aws/SAMN001.fa.gz"),
        ("SAMN002", "SRR002", "GCA_002", 1, "Escherichia coli", "PASS",
         "batch_001.tar.xz", "https://osf.io/batch_001", "https://aws/SAMN002.fa.gz"),
        ("SAMN003", "SRR003", "GCA_003", 1, "Staphylococcus aureus", "PASS",
         "batch_002.tar.xz", "https://osf.io/batch_002", "https://aws/SAMN003.fa.gz"),
        ("SAMN004", "SRR004", "GCA_004", 1, "Klebsiella pneumoniae", "PASS",
         "batch_002.tar.xz", "https://osf.io/batch_002", "https://aws/SAMN004.fa.gz"),
    ]
    conn.executemany("INSERT INTO assembly VALUES (?,?,?,?,?,?,?,?,?)", assemblies)

    runs = [
        ("SRR001", "SAMN001", 1),
        ("SRR002", "SAMN002", 1),
        ("SRR003", "SAMN003", 1),
        ("SRR004", "SAMN004", 1),
    ]
    conn.executemany("INSERT INTO run VALUES (?,?,?)", runs)

    ena = [
        ("SRR001", "SAMN001", "United Kingdom", "2022-03-15", "Homo sapiens", "blood"),
        ("SRR002", "SAMN002", "Germany", "2023-06-01", "Homo sapiens", "stool"),
        ("SRR003", "SAMN003", "United Kingdom", "2023-01-10", "Homo sapiens", "blood"),
        ("SRR004", "SAMN004", "Japan", "2020-08-05", "Homo sapiens", "wound"),
    ]
    conn.executemany(
        "INSERT INTO ena_202505_used "
        "(run_accession, sample_accession, country, collection_date, host, isolation_source) "
        "VALUES (?,?,?,?,?,?)",
        ena,
    )

    checkm2 = [
        ("SAMN001", 99.5, 0.5, 5100000, 50.5, 150000),
        ("SAMN002", 97.2, 1.2, 4900000, 50.8, 120000),
        ("SAMN003", 95.5, 2.1, 2800000, 32.8, 80000),
        ("SAMN004", 99.1, 0.3, 5600000, 57.3, 200000),
    ]
    conn.executemany("INSERT INTO checkm2 VALUES (?,?,?,?,?,?)", checkm2)

    conn.commit()
    conn.close()
    return db_path


# ---------------------------------------------------------------------------
# Feature 1: Species typo suggestions
# ---------------------------------------------------------------------------


class TestSuggestSpecies:
    def test_finds_close_match(self):
        available = ["Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae"]
        suggestions = _suggest_species("Escherichia coli", available)
        assert "Escherichia coli" in suggestions

    def test_suggests_typo(self):
        available = ["Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae"]
        suggestions = _suggest_species("Eschericia coli", available)  # typo: missing h
        assert len(suggestions) > 0
        assert "Escherichia coli" in suggestions

    def test_no_match_returns_empty(self):
        available = ["Escherichia coli", "Staphylococcus aureus"]
        suggestions = _suggest_species("Completely different organism xyz", available)
        assert suggestions == []

    def test_cli_shows_suggestions_for_bad_species(self, runner, test_db):
        """When --species has no results, CLI should suggest alternatives."""
        result = runner.invoke(
            main,
            [
                "query",
                "--db-path", str(test_db),
                "--species", "Eschericia coli",  # typo
            ],
        )
        # Should fail with exit code 1
        assert result.exit_code == 1
        # Should suggest correct spelling
        assert "Escherichia coli" in result.output or "Did you mean" in result.output


# ---------------------------------------------------------------------------
# Feature 2: Terminal-aware output formatting
# ---------------------------------------------------------------------------


class TestOutputFormatting:
    def test_resolve_format_auto_non_tty(self, monkeypatch):
        """In non-TTY context, auto should resolve to tsv."""
        monkeypatch.setattr("sys.stdout.isatty", lambda: False)
        assert resolve_format("auto") == "tsv"

    def test_resolve_format_auto_tty(self, monkeypatch):
        """In TTY context, auto should resolve to table."""
        monkeypatch.setattr("sys.stdout.isatty", lambda: True)
        assert resolve_format("auto") == "table"

    def test_resolve_format_explicit(self):
        for fmt in ("table", "tsv", "csv", "json"):
            assert resolve_format(fmt) == fmt

    def test_print_tsv(self, capsys):
        df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
        print_dataframe(df, "tsv")
        captured = capsys.readouterr()
        assert "a\tb\n" in captured.out
        assert "1\tx" in captured.out

    def test_print_csv(self, capsys):
        df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
        print_dataframe(df, "csv")
        captured = capsys.readouterr()
        assert "a,b" in captured.out

    def test_print_json(self, capsys):
        df = pd.DataFrame({"a": [1], "b": ["x"]})
        print_dataframe(df, "json")
        captured = capsys.readouterr()
        assert '"a"' in captured.out
        assert '"b"' in captured.out

    def test_print_table(self, capsys):
        df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
        print_dataframe(df, "table")
        captured = capsys.readouterr()
        assert "a" in captured.out and "b" in captured.out

    def test_format_option_in_query_cmd_help(self, runner):
        result = runner.invoke(main, ["query", "--help"])
        assert "--format" in result.output

    def test_format_option_in_list_species_help(self, runner):
        result = runner.invoke(main, ["list-species", "--help"])
        assert "--format" in result.output

    def test_format_option_in_mlst_query_help(self, runner):
        result = runner.invoke(main, ["mlst-query", "--help"])
        assert "--format" in result.output


# ---------------------------------------------------------------------------
# Feature 3: info subcommand
# ---------------------------------------------------------------------------


class TestInfoSubcommand:
    def test_info_help(self, runner):
        result = runner.invoke(main, ["info", "--help"])
        assert result.exit_code == 0
        assert "SAMPLE_ACCESSION" in result.output

    def test_info_known_sample(self, runner, test_db):
        result = runner.invoke(
            main,
            ["info", "SAMN001", "--db-path", str(test_db)],
        )
        assert result.exit_code == 0
        assert "SAMN001" in result.output
        assert "Escherichia coli" in result.output
        assert "Assembly" in result.output

    def test_info_shows_checkm2(self, runner, test_db):
        result = runner.invoke(
            main,
            ["info", "SAMN001", "--db-path", str(test_db)],
        )
        assert "Completeness" in result.output
        assert "Contamination" in result.output

    def test_info_shows_ena(self, runner, test_db):
        result = runner.invoke(
            main,
            ["info", "SAMN001", "--db-path", str(test_db)],
        )
        assert result.exit_code == 0
        assert "ENA" in result.output
        assert "United Kingdom" in result.output

    def test_info_unknown_sample_exits_nonzero(self, runner, test_db):
        result = runner.invoke(
            main,
            ["info", "NONEXISTENT", "--db-path", str(test_db)],
        )
        assert result.exit_code != 0

    def test_info_no_db_exits_nonzero(self, runner, tmp_path):
        result = runner.invoke(
            main,
            ["info", "SAMN001", "--cache-dir", str(tmp_path)],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Feature 4: Dry-run
# ---------------------------------------------------------------------------


class TestDryRun:
    def test_accessions_dry_run(self, runner, tmp_path):
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("SAMN001\nSAMN002\n")
        result = runner.invoke(
            main,
            [
                "accessions",
                str(acc_file),
                "--output", str(tmp_path / "out"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry-run" in result.output
        assert "SAMN001" in result.output
        assert "SAMN002" in result.output
        # Output dir should not be created
        assert not (tmp_path / "out").exists()

    def test_query_dry_run(self, runner, test_db, tmp_path):
        result = runner.invoke(
            main,
            [
                "query",
                "--db-path", str(test_db),
                "--species", "Escherichia coli",
                "--output", str(tmp_path / "out"),
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert "dry-run" in result.output
        # Nothing should have been downloaded
        assert not (tmp_path / "out").exists() or not any((tmp_path / "out").iterdir())

    def test_query_dry_run_help(self, runner):
        result = runner.invoke(main, ["query", "--help"])
        assert "--dry-run" in result.output

    def test_accessions_dry_run_help(self, runner):
        result = runner.invoke(main, ["accessions", "--help"])
        assert "--dry-run" in result.output


# ---------------------------------------------------------------------------
# Feature 5: mlst-query subcommand
# ---------------------------------------------------------------------------


class TestMlstQuerySubcommand:
    def test_mlst_query_help(self, runner):
        result = runner.invoke(main, ["mlst-query", "--help"])
        assert result.exit_code == 0
        assert "--st" in result.output
        assert "--scheme" in result.output
        assert "--species" in result.output

    def test_mlst_query_subcommand_exists(self, runner):
        result = runner.invoke(main, ["--help"])
        assert "mlst-query" in result.output


# ---------------------------------------------------------------------------
# Feature 6: TOML filter files
# ---------------------------------------------------------------------------


class TestTomlFilterFiles:
    def test_load_toml_filters_none(self, tmp_path):
        result = _load_toml_filters(None)
        assert result == {}

    def test_load_toml_filters_missing_file(self, tmp_path):
        result = _load_toml_filters(tmp_path / "nonexistent.toml")
        assert result == {}

    def test_load_toml_filters_valid(self, tmp_path):
        toml_file = tmp_path / "filters.toml"
        toml_file.write_text(
            '[filters]\nspecies = "Escherichia coli"\ncountry = "United Kingdom"\n'
        )
        result = _load_toml_filters(toml_file)
        assert result["species"] == "Escherichia coli"
        assert result["country"] == "United Kingdom"

    def test_load_toml_filters_with_numbers(self, tmp_path):
        toml_file = tmp_path / "filters.toml"
        toml_file.write_text("[query]\nyear_from = 2020\nyear_to = 2023\nn = 500\n")
        result = _load_toml_filters(toml_file)
        assert result["year_from"] == 2020
        assert result["year_to"] == 2023

    def test_query_uses_toml_filter(self, runner, test_db, tmp_path):
        toml_file = tmp_path / "filters.toml"
        toml_file.write_text('[filters]\nspecies = "Escherichia coli"\n')
        result = runner.invoke(
            main,
            [
                "query",
                "--db-path", str(test_db),
                "--filter", str(toml_file),
            ],
        )
        assert result.exit_code == 0
        assert "SAMN001" in result.output or "SAMN002" in result.output

    def test_cli_flag_overrides_toml(self, runner, test_db, tmp_path):
        """CLI --species should override TOML species."""
        toml_file = tmp_path / "filters.toml"
        toml_file.write_text('[filters]\nspecies = "Nonexistent species 999"\n')
        result = runner.invoke(
            main,
            [
                "query",
                "--db-path", str(test_db),
                "--filter", str(toml_file),
                "--species", "Escherichia coli",  # should override TOML
            ],
        )
        assert result.exit_code == 0
        assert "SAMN001" in result.output or "SAMN002" in result.output


# ---------------------------------------------------------------------------
# Feature 7: summarise subcommand
# ---------------------------------------------------------------------------


class TestSummariseSubcommand:
    def test_summarise_help(self, runner):
        result = runner.invoke(main, ["summarise", "--help"])
        assert result.exit_code == 0
        assert "--by" in result.output
        assert "--top" in result.output

    def test_summarise_from_tsv(self, runner, tmp_path):
        tsv_file = tmp_path / "data.tsv"
        df = pd.DataFrame({
            "sample": ["SAMN001", "SAMN002", "SAMN003"],
            "species": ["Escherichia coli", "Escherichia coli", "Staphylococcus aureus"],
            "country": ["UK", "UK", "Japan"],
        })
        df.to_csv(tsv_file, sep="\t", index=False)
        result = runner.invoke(
            main,
            ["summarise", "--by", "species", "--input", str(tsv_file)],
        )
        assert result.exit_code == 0
        assert "Escherichia coli" in result.output
        assert "Staphylococcus aureus" in result.output

    def test_summarise_top_n(self, runner, tmp_path):
        tsv_file = tmp_path / "data.tsv"
        df = pd.DataFrame({
            "sample": [f"S{i}" for i in range(10)],
            "species": ["A"] * 5 + ["B"] * 3 + ["C"] * 2,
        })
        df.to_csv(tsv_file, sep="\t", index=False)
        result = runner.invoke(
            main,
            ["summarise", "--by", "species", "--top", "2", "--input", str(tsv_file)],
        )
        assert result.exit_code == 0
        # Should have A and B but not C
        assert "A" in result.output
        assert "B" in result.output
        assert "C" not in result.output

    def test_summarise_by_country(self, runner, tmp_path):
        tsv_file = tmp_path / "data.tsv"
        df = pd.DataFrame({
            "sample": ["SAMN001", "SAMN002", "SAMN003"],
            "species": ["Escherichia coli"] * 3,
            "country": ["UK", "UK", "Japan"],
        })
        df.to_csv(tsv_file, sep="\t", index=False)
        result = runner.invoke(
            main,
            ["summarise", "--by", "country", "--input", str(tsv_file)],
        )
        assert result.exit_code == 0
        assert "UK" in result.output
        assert "Japan" in result.output

    def test_summarise_bad_column(self, runner, tmp_path):
        tsv_file = tmp_path / "data.tsv"
        df = pd.DataFrame({"sample": ["SAMN001"], "species": ["Ecoli"]})
        df.to_csv(tsv_file, sep="\t", index=False)
        result = runner.invoke(
            main,
            ["summarise", "--by", "nonexistent_col", "--input", str(tsv_file)],
        )
        assert result.exit_code != 0

    def test_summarise_with_db(self, runner, test_db):
        result = runner.invoke(
            main,
            ["summarise", "--by", "species", "--db-path", str(test_db)],
        )
        assert result.exit_code == 0
        assert "Escherichia coli" in result.output


# ---------------------------------------------------------------------------
# Feature 8: Config file
# ---------------------------------------------------------------------------


class TestConfigSubcommand:
    def test_config_help(self, runner):
        result = runner.invoke(main, ["config", "--help"])
        assert result.exit_code == 0
        assert "init" in result.output
        assert "show" in result.output
        assert "set" in result.output

    def test_config_init(self, runner, tmp_path):
        cfg_path = tmp_path / "config.toml"
        # Write config to a tmp path by monkeypatching CONFIG_PATH
        import atbfetcher.config as cfg_mod
        orig = cfg_mod.CONFIG_PATH
        cfg_mod.CONFIG_PATH = cfg_path
        try:
            result = runner.invoke(main, ["config", "init"])
            assert result.exit_code == 0
            assert cfg_path.exists()
        finally:
            cfg_mod.CONFIG_PATH = orig

    def test_config_show_no_file(self, runner, tmp_path):
        import atbfetcher.config as cfg_mod
        orig = cfg_mod.CONFIG_PATH
        cfg_mod.CONFIG_PATH = tmp_path / "config.toml"
        try:
            result = runner.invoke(main, ["config", "show"])
            assert result.exit_code == 0
            assert "No config file" in result.output
        finally:
            cfg_mod.CONFIG_PATH = orig

    def test_write_and_load_config(self, tmp_path):
        cfg_path = tmp_path / "config.toml"
        cfg = {"defaults": {"output_format": "tsv", "cache_dir": str(tmp_path)}}
        write_config(cfg, path=cfg_path)
        loaded = load_config(path=cfg_path)
        assert loaded["defaults"]["output_format"] == "tsv"
        assert loaded["defaults"]["cache_dir"] == str(tmp_path)

    def test_default_config_structure(self):
        cfg = default_config()
        assert "defaults" in cfg
        assert "output_format" in cfg["defaults"]
        assert "cache_dir" in cfg["defaults"]

    def test_config_set_and_get(self, tmp_path):
        cfg_path = tmp_path / "config.toml"
        write_config(default_config(), path=cfg_path)

        import atbfetcher.config as cfg_mod
        orig = cfg_mod.CONFIG_PATH
        cfg_mod.CONFIG_PATH = cfg_path
        try:
            runner = CliRunner()
            result = runner.invoke(main, ["config", "set", "output_format", "json"])
            assert result.exit_code == 0

            loaded = load_config(path=cfg_path)
            assert loaded["defaults"]["output_format"] == "json"
        finally:
            cfg_mod.CONFIG_PATH = orig

    def test_config_set_numeric(self, tmp_path):
        cfg_path = tmp_path / "config.toml"
        write_config(default_config(), path=cfg_path)

        import atbfetcher.config as cfg_mod
        orig = cfg_mod.CONFIG_PATH
        cfg_mod.CONFIG_PATH = cfg_path
        try:
            runner = CliRunner()
            result = runner.invoke(main, ["config", "set", "threads", "8"])
            assert result.exit_code == 0

            loaded = load_config(path=cfg_path)
            assert loaded["defaults"]["threads"] == 8
        finally:
            cfg_mod.CONFIG_PATH = orig


# ---------------------------------------------------------------------------
# Regression: existing commands still work after changes
# ---------------------------------------------------------------------------


class TestExistingCommandsUnchanged:
    def test_species_help(self, runner):
        result = runner.invoke(main, ["species", "--help"])
        assert result.exit_code == 0
        assert "SPECIES_NAME" in result.output

    def test_mlst_help(self, runner):
        result = runner.invoke(main, ["mlst", "--help"])
        assert result.exit_code == 0

    def test_accessions_help(self, runner):
        result = runner.invoke(main, ["accessions", "--help"])
        assert result.exit_code == 0

    def test_query_help(self, runner):
        result = runner.invoke(main, ["query", "--help"])
        assert result.exit_code == 0

    def test_list_species_help(self, runner):
        result = runner.invoke(main, ["list-species", "--help"])
        assert result.exit_code == 0

    def test_species_count_help(self, runner):
        result = runner.invoke(main, ["species-count", "--help"])
        assert result.exit_code == 0

    def test_list_countries_help(self, runner):
        result = runner.invoke(main, ["list-countries", "--help"])
        assert result.exit_code == 0

    def test_list_hosts_help(self, runner):
        result = runner.invoke(main, ["list-hosts", "--help"])
        assert result.exit_code == 0

    def test_download_db_help(self, runner):
        result = runner.invoke(main, ["download-db", "--help"])
        assert result.exit_code == 0
