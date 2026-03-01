"""Tests for species name handling."""

import pandas as pd

from atbfetcher.species import (
    clean_species_name,
    get_samples_for_species,
    is_placeholder_species,
    list_species,
    normalize_for_matching,
)


class TestIsPlaceholderSpecies:
    """Tests for GTDB placeholder species detection."""

    def test_detects_numeric_placeholder(self):
        assert is_placeholder_species("Bacillus sp000746275") is True

    def test_detects_various_placeholders(self):
        assert is_placeholder_species("Clostridium sp001234567") is True
        assert is_placeholder_species("Enterobacter sp900112345") is True

    def test_real_species_not_placeholder(self):
        assert is_placeholder_species("Escherichia coli") is False
        assert is_placeholder_species("Staphylococcus aureus") is False

    def test_species_with_suffix_not_placeholder(self):
        assert is_placeholder_species("Enterobacter hormaechei_A") is False

    def test_genus_only_not_placeholder(self):
        assert is_placeholder_species("Pseudomonas") is False


class TestCleanSpeciesName:
    """Tests for GTDB suffix stripping."""

    def test_strips_single_letter_suffix(self):
        assert clean_species_name("Enterobacter hormaechei_A") == "Enterobacter hormaechei"

    def test_strips_various_suffixes(self):
        assert clean_species_name("Klebsiella variicola_B") == "Klebsiella variicola"
        assert clean_species_name("Citrobacter freundii_C") == "Citrobacter freundii"

    def test_strips_genus_level_suffix(self):
        assert clean_species_name("Campylobacter_D jejuni") == "Campylobacter jejuni"
        assert clean_species_name("Clostridium_B difficile") == "Clostridium difficile"

    def test_no_suffix_unchanged(self):
        assert clean_species_name("Escherichia coli") == "Escherichia coli"

    def test_single_word_genus_only(self):
        assert clean_species_name("Pseudomonas") == "Pseudomonas"

    def test_does_not_strip_lowercase_suffix(self):
        # Only uppercase single-letter suffixes should be stripped
        assert clean_species_name("Foo bar_a") == "Foo bar_a"

    def test_does_not_strip_multi_letter_suffix(self):
        assert clean_species_name("Foo bar_AB") == "Foo bar_AB"


class TestNormalizeForMatching:
    """Tests for cross-dataset name normalization."""

    def test_lowercases_and_strips(self):
        assert normalize_for_matching("Escherichia coli") == "escherichia coli"

    def test_strips_suffix_and_lowercases(self):
        assert normalize_for_matching("Enterobacter hormaechei_A") == "enterobacter hormaechei"


class TestListSpecies:
    """Tests for listing unique species."""

    def test_unique_sorted(self, sample_species_calls_df):
        names = list_species(sample_species_calls_df)
        # "Enterobacter hormaechei_A" cleaned to "Enterobacter hormaechei"
        assert "Enterobacter hormaechei" in names
        assert "Escherichia coli" in names
        assert names == sorted(names)

    def test_no_duplicates_after_cleaning(self, sample_species_calls_df):
        names = list_species(sample_species_calls_df)
        assert len(names) == len(set(names))

    def test_excludes_placeholder_species(self, sample_species_calls_df):
        names = list_species(sample_species_calls_df)
        assert not any("sp000" in n for n in names)
        assert "Bacillus sp000746275" not in names

    def test_raw_mode_keeps_suffixes(self, sample_species_calls_df):
        names = list_species(sample_species_calls_df, raw=True)
        assert "Enterobacter hormaechei_A" in names

    def test_raw_mode_excludes_placeholders(self, sample_species_calls_df):
        names = list_species(sample_species_calls_df, raw=True)
        assert "Bacillus sp000746275" not in names


class TestGetSamplesForSpecies:
    """Tests for filtering samples by species."""

    def test_returns_matching_samples(self, sample_species_calls_df):
        result = get_samples_for_species("Escherichia coli", sample_species_calls_df)
        assert len(result) == 3
        assert set(result["sample"]) == {"SAMN001", "SAMN002", "SAMN007"}

    def test_matches_cleaned_name(self, sample_species_calls_df):
        # Searching for "Enterobacter hormaechei" should match "_A" suffixed entries
        result = get_samples_for_species(
            "Enterobacter hormaechei", sample_species_calls_df
        )
        assert len(result) == 2
        assert set(result["sample"]) == {"SAMN003", "SAMN004"}

    def test_no_match_returns_empty(self, sample_species_calls_df):
        result = get_samples_for_species("Nonexistent species", sample_species_calls_df)
        assert len(result) == 0

    def test_case_insensitive_match(self, sample_species_calls_df):
        result = get_samples_for_species("escherichia coli", sample_species_calls_df)
        assert len(result) == 3

    def test_placeholder_species_returns_empty(self, sample_species_calls_df):
        result = get_samples_for_species("Bacillus sp000746275", sample_species_calls_df)
        assert len(result) == 0
