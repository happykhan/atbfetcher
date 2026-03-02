"""Tests for SQLite metadata query module."""

import sqlite3

import pytest

from atbfetcher.query import find_sqlite_db, list_countries, list_hosts, query_metadata


@pytest.fixture
def test_db(tmp_path):
    """Create a small test SQLite database mimicking ATB metadata."""
    db_path = tmp_path / "test.sqlite"
    conn = sqlite3.connect(str(db_path))

    # Create tables matching the ATB schema
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

    # Insert test data - assembly
    assemblies = [
        (
            "SAMN001",
            "SRR001",
            "GCA_001",
            1,
            "Escherichia coli",
            "PASS",
            "batch_001.tar.xz",
            "https://osf.io/batch_001",
            "https://aws/SAMN001.fa.gz",
        ),
        (
            "SAMN002",
            "SRR002",
            "GCA_002",
            1,
            "Escherichia coli",
            "PASS",
            "batch_001.tar.xz",
            "https://osf.io/batch_001",
            "https://aws/SAMN002.fa.gz",
        ),
        (
            "SAMN003",
            "SRR003",
            "GCA_003",
            1,
            "Escherichia coli",
            "PASS",
            "batch_001.tar.xz",
            "https://osf.io/batch_001",
            "https://aws/SAMN003.fa.gz",
        ),
        (
            "SAMN004",
            "SRR004",
            "GCA_004",
            1,
            "Staphylococcus aureus",
            "PASS",
            "batch_002.tar.xz",
            "https://osf.io/batch_002",
            "https://aws/SAMN004.fa.gz",
        ),
        (
            "SAMN005",
            "SRR005",
            "GCA_005",
            1,
            "Staphylococcus aureus",
            "PASS",
            "batch_002.tar.xz",
            "https://osf.io/batch_002",
            "https://aws/SAMN005.fa.gz",
        ),
        (
            "SAMN006",
            "SRR006",
            "GCA_006",
            1,
            "Escherichia coli",
            "CHECKM2_MIN_COMPL",
            "batch_003.tar.xz",
            "https://osf.io/batch_003",
            "https://aws/SAMN006.fa.gz",
        ),
        (
            "SAMN007",
            "SRR007",
            "GCA_007",
            0,
            "Klebsiella pneumoniae",
            "PASS",
            None,
            None,
            "https://aws/SAMN007.fa.gz",
        ),
        (
            "SAMN008",
            "SRR008",
            "GCA_008",
            1,
            "Escherichia coli_A",
            "PASS",
            "batch_001.tar.xz",
            "https://osf.io/batch_001",
            "https://aws/SAMN008.fa.gz",
        ),
    ]
    conn.executemany("INSERT INTO assembly VALUES (?,?,?,?,?,?,?,?,?)", assemblies)

    # Insert test data - run
    runs = [
        ("SRR001", "SAMN001", 1),
        ("SRR002", "SAMN002", 1),
        ("SRR003", "SAMN003", 1),
        ("SRR004", "SAMN004", 1),
        ("SRR005", "SAMN005", 1),
        ("SRR006", "SAMN006", 1),
        ("SRR007", "SAMN007", 1),
        ("SRR008", "SAMN008", 1),
    ]
    conn.executemany("INSERT INTO run VALUES (?,?,?)", runs)

    # Insert test data - ENA metadata
    ena = [
        ("SRR001", "SAMN001", "United Kingdom", "2022-03-15", "Homo sapiens", "blood"),
        ("SRR002", "SAMN002", "United Kingdom:England", "2023-06-01", "Homo sapiens", "stool"),
        ("SRR003", "SAMN003", "Germany", "2021-11-20", "Homo sapiens", "urine"),
        ("SRR004", "SAMN004", "United Kingdom", "2023-01-10", "Homo sapiens", "blood"),
        ("SRR005", "SAMN005", "Japan", "2020-08-05", "Homo sapiens", "wound"),
        ("SRR006", "SAMN006", "United Kingdom", "2022-07-30", "Homo sapiens", "blood"),
        ("SRR007", "SAMN007", "Germany", "2023-04-12", "Gallus gallus", "caecum"),
        ("SRR008", "SAMN008", "France", "2023-09-01", "Homo sapiens", "stool"),
    ]
    conn.executemany(
        "INSERT INTO ena_202505_used (run_accession, sample_accession, country, "
        "collection_date, host, isolation_source) VALUES (?,?,?,?,?,?)",
        ena,
    )

    # Insert test data - CheckM2
    checkm2 = [
        ("SAMN001", 99.5, 0.5, 5100000, 50.5, 150000),
        ("SAMN002", 97.2, 1.2, 4900000, 50.8, 120000),
        ("SAMN003", 98.0, 0.8, 5200000, 50.2, 95000),
        ("SAMN004", 95.5, 2.1, 2800000, 32.8, 80000),
        ("SAMN005", 93.0, 3.5, 2900000, 33.1, 60000),
        ("SAMN006", 60.0, 15.0, 5300000, 51.0, 10000),
        ("SAMN007", 99.1, 0.3, 5600000, 57.3, 200000),
        ("SAMN008", 96.0, 1.0, 5050000, 50.2, 130000),
    ]
    conn.executemany("INSERT INTO checkm2 VALUES (?,?,?,?,?,?)", checkm2)

    conn.commit()
    conn.close()
    return db_path


class TestFindSqliteDb:
    """Tests for find_sqlite_db."""

    def test_finds_explicit_path(self, test_db, tmp_path):
        result = find_sqlite_db(test_db, tmp_path)
        assert result == test_db

    def test_finds_in_cache_dir(self, test_db, tmp_path):
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        db_in_cache = cache_dir / "atb.metadata.202505.sqlite"
        db_in_cache.symlink_to(test_db)
        result = find_sqlite_db(None, cache_dir)
        assert result == db_in_cache

    def test_returns_none_when_not_found(self, tmp_path):
        result = find_sqlite_db(None, tmp_path / "nonexistent")
        assert result is None


class TestQueryMetadata:
    """Tests for query_metadata."""

    def test_query_all_hq(self, test_db):
        """Query all HQ genomes returns only PASS samples on OSF."""
        results = query_metadata(test_db, hq_only=True)
        assert len(results) == 6  # SAMN001-005 + SAMN008 (all PASS + on OSF)
        assert "SAMN006" not in results["sample"].values  # not HQ
        assert "SAMN007" not in results["sample"].values  # not on OSF

    def test_query_no_hq_filter(self, test_db):
        """Without HQ filter, includes non-PASS genomes."""
        results = query_metadata(test_db, hq_only=False)
        assert "SAMN006" in results["sample"].values  # not HQ but included
        assert "SAMN007" not in results["sample"].values  # still excluded: not on OSF

    def test_filter_by_species(self, test_db):
        """Filter by exact species name."""
        results = query_metadata(test_db, species="Escherichia coli")
        # Should match "Escherichia coli" and "Escherichia coli_A" (GTDB suffix)
        species_set = set(results["species"])
        assert "Escherichia coli" in species_set
        assert "Escherichia coli_A" in species_set
        # Only PASS + on OSF E. coli: SAMN001, SAMN002, SAMN003, SAMN008
        assert len(results) == 4

    def test_filter_by_species_staph(self, test_db):
        """Filter S. aureus returns only S. aureus genomes."""
        results = query_metadata(test_db, species="Staphylococcus aureus")
        assert len(results) == 2
        assert set(results["sample"]) == {"SAMN004", "SAMN005"}

    def test_filter_by_country(self, test_db):
        """Country filter uses prefix matching."""
        results = query_metadata(test_db, country="United Kingdom")
        # Should match "United Kingdom" and "United Kingdom:England"
        assert len(results) >= 2
        for _, row in results.iterrows():
            assert row["country"].startswith("United Kingdom")

    def test_filter_by_country_and_species(self, test_db):
        """Combine species and country filters."""
        results = query_metadata(test_db, species="Escherichia coli", country="United Kingdom")
        for _, row in results.iterrows():
            assert row["species"].startswith("Escherichia coli")
            assert row["country"].startswith("United Kingdom")

    def test_filter_by_year_range(self, test_db):
        """Year range filter works on collection_date."""
        results = query_metadata(test_db, year_from=2023, year_to=2023)
        for _, row in results.iterrows():
            assert row["collection_date"].startswith("2023")

    def test_filter_by_year_from(self, test_db):
        """Year-from filter (no upper bound)."""
        results = query_metadata(test_db, year_from=2023)
        for _, row in results.iterrows():
            year = int(row["collection_date"][:4])
            assert year >= 2023

    def test_filter_by_host(self, test_db):
        """Host filter with LIKE match."""
        results = query_metadata(test_db, host="Gallus")
        assert len(results) == 0  # SAMN007 has Gallus but is not on OSF

    def test_filter_by_isolation_source(self, test_db):
        """Isolation source filter with LIKE match."""
        results = query_metadata(test_db, isolation_source="blood")
        assert len(results) >= 1
        for _, row in results.iterrows():
            assert "blood" in row["isolation_source"].lower()

    def test_limit_results(self, test_db):
        """Limit number of results."""
        results = query_metadata(test_db, limit=2)
        assert len(results) == 2

    def test_checkm2_min_completeness(self, test_db):
        """CheckM2 completeness filter."""
        results = query_metadata(test_db, hq_only=False, min_completeness=95.0)
        for _, row in results.iterrows():
            assert row["Completeness_Specific"] >= 95.0

    def test_checkm2_max_contamination(self, test_db):
        """CheckM2 contamination filter."""
        results = query_metadata(test_db, hq_only=False, max_contamination=2.0)
        for _, row in results.iterrows():
            assert row["Contamination"] <= 2.0

    def test_checkm2_genome_size_range(self, test_db):
        """CheckM2 genome size range filter."""
        results = query_metadata(test_db, min_genome_size=4_000_000, max_genome_size=5_500_000)
        for _, row in results.iterrows():
            assert 4_000_000 <= row["Genome_Size"] <= 5_500_000

    def test_combined_filters(self, test_db):
        """Multiple filters combined."""
        results = query_metadata(
            test_db,
            species="Escherichia coli",
            country="United Kingdom",
            year_from=2022,
            year_to=2023,
            isolation_source="blood",
        )
        # SAMN001: E. coli, UK, 2022, blood -> match
        assert len(results) == 1
        assert results.iloc[0]["sample"] == "SAMN001"

    def test_empty_result(self, test_db):
        """Query with impossible filters returns empty."""
        results = query_metadata(test_db, species="Nonexistent species")
        assert len(results) == 0

    def test_columns_with_ena(self, test_db):
        """Results include ENA columns when ENA filters are used."""
        results = query_metadata(test_db, country="Germany")
        assert "country" in results.columns
        assert "collection_date" in results.columns
        assert "host" in results.columns

    def test_columns_without_ena(self, test_db):
        """Results don't include ENA columns when no ENA filters are used."""
        results = query_metadata(test_db, species="Escherichia coli")
        assert "country" not in results.columns

    def test_columns_with_checkm2(self, test_db):
        """Results include CheckM2 columns when CheckM2 filters are used."""
        results = query_metadata(test_db, min_completeness=90.0)
        assert "Completeness_Specific" in results.columns
        assert "Contamination" in results.columns
        assert "Genome_Size" in results.columns


class TestListCountries:
    """Tests for list_countries."""

    def test_lists_all_countries(self, test_db):
        countries = list_countries(test_db)
        assert "Germany" in countries
        assert "Japan" in countries
        assert "United Kingdom" in countries

    def test_lists_countries_for_species(self, test_db):
        countries = list_countries(test_db, species="Staphylococcus aureus")
        assert "United Kingdom" in countries
        assert "Japan" in countries
        # E. coli-only countries should not appear
        assert "Germany" not in countries  # Only SAMN003 (E. coli) is in Germany


class TestListHosts:
    """Tests for list_hosts."""

    def test_lists_all_hosts(self, test_db):
        hosts = list_hosts(test_db)
        assert "Homo sapiens" in hosts

    def test_lists_hosts_for_species(self, test_db):
        hosts = list_hosts(test_db, species="Escherichia coli")
        assert "Homo sapiens" in hosts
