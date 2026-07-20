#!/usr/bin/env python3

import os
import shutil
import tempfile

import pytest

from kegg_pathways_completeness.bin.utils import (
    load_modules_kos,
    parse_module_kos,
    sanity_check_give_completeness_output,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs"""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def modules_table_file(temp_dir):
    path = os.path.join(temp_dir, "modules_table.tsv")
    with open(path, "w") as f:
        f.write("module\tdefinition\tname\tclass\n")
        f.write("M00001\t(K00844,K12407) K01810\tGlycolysis\tCarbohydrate metabolism\n")
        f.write("M00002\tK01647 K01681\tTCA cycle\tCarbohydrate metabolism\n")
    return path


def write_output_table(temp_dir, rows, header=None):
    path = os.path.join(temp_dir, "output.tsv")
    header = header or [
        "module_accession",
        "completeness",
        "pathway_name",
        "pathway_class",
        "matching_ko",
        "missing_ko",
    ]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    return path


class TestParseModuleKos:
    def test_extracts_kos_from_simple_definition(self):
        assert parse_module_kos("K00134 K00927") == {"K00134", "K00927"}

    def test_extracts_kos_with_parentheses_and_operators(self):
        definition = "(K00844,K12407) (((K00134,K00150) K00927),K11389)"
        assert parse_module_kos(definition) == {
            "K00844",
            "K12407",
            "K00134",
            "K00150",
            "K00927",
            "K11389",
        }

    def test_ignores_weight_suffix(self):
        # matching_ko/missing_ko with --include-weights look like "K00134(0.5)"
        assert parse_module_kos("K00134(0.5),K00927(1.0)") == {"K00134", "K00927"}

    def test_no_kos_found(self):
        assert parse_module_kos("") == set()


class TestLoadModulesKos:
    def test_loads_kos_per_module(self, modules_table_file):
        modules_kos = load_modules_kos(modules_table_file)
        assert modules_kos["M00001"] == {"K00844", "K12407", "K01810"}
        assert modules_kos["M00002"] == {"K01647", "K01681"}

    def test_missing_required_column_raises(self, temp_dir):
        path = os.path.join(temp_dir, "bad_modules_table.tsv")
        with open(path, "w") as f:
            f.write("module\tname\tclass\n")
            f.write("M00001\tGlycolysis\tCarbohydrate metabolism\n")

        with pytest.raises(ValueError):
            load_modules_kos(path)


class TestSanityCheckGiveCompletenessOutput:
    def test_consistent_output_has_no_errors(self, temp_dir, modules_table_file):
        output_file = write_output_table(
            temp_dir,
            rows=[
                [
                    "M00001",
                    "66.67",
                    "Glycolysis",
                    "Carbohydrate metabolism",
                    "K00844,K12407",
                    "K01810",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert errors == []

    def test_flags_matching_ko_not_in_definition(self, temp_dir, modules_table_file):
        output_file = write_output_table(
            temp_dir,
            rows=[
                [
                    "M00001",
                    "100.0",
                    "Glycolysis",
                    "Carbohydrate metabolism",
                    "K00844,K99999",
                    "",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert len(errors) == 1
        assert "matching_ko K99999" in errors[0]
        assert "M00001" in errors[0]

    def test_flags_missing_ko_not_in_definition(self, temp_dir, modules_table_file):
        output_file = write_output_table(
            temp_dir,
            rows=[
                [
                    "M00002",
                    "50.0",
                    "TCA cycle",
                    "Carbohydrate metabolism",
                    "K01647",
                    "K99999",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert len(errors) == 1
        assert "missing_ko K99999" in errors[0]
        assert "M00002" in errors[0]

    def test_flags_unknown_module(self, temp_dir, modules_table_file):
        output_file = write_output_table(
            temp_dir,
            rows=[
                [
                    "M09999",
                    "100.0",
                    "Unknown",
                    "Unknown",
                    "K00844",
                    "",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert len(errors) == 1
        assert "M09999" in errors[0]
        assert "not found" in errors[0]

    def test_handles_per_contig_output_with_contig_column(
        self, temp_dir, modules_table_file
    ):
        output_file = write_output_table(
            temp_dir,
            header=[
                "contig",
                "module_accession",
                "completeness",
                "pathway_name",
                "pathway_class",
                "matching_ko",
                "missing_ko",
            ],
            rows=[
                [
                    "contig1",
                    "M00001",
                    "66.67",
                    "Glycolysis",
                    "Carbohydrate metabolism",
                    "K00844,K12407",
                    "K01810",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert errors == []

    def test_handles_weighted_ko_format(self, temp_dir, modules_table_file):
        output_file = write_output_table(
            temp_dir,
            rows=[
                [
                    "M00001",
                    "66.67",
                    "Glycolysis",
                    "Carbohydrate metabolism",
                    "K00844(0.5),K12407(0.5)",
                    "K01810(1.0)",
                ],
            ],
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert errors == []


class TestSanityCheckOnRealFixtures:
    """Cross-check against the real give_completeness fixtures/outputs."""

    @pytest.mark.parametrize(
        "output_filename",
        [
            "test_combined_pathways.tsv",
            "test_combined_contigs.tsv",
        ],
    )
    def test_real_output_is_consistent(self, output_filename):
        repo_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        output_file = os.path.join(
            repo_root, "tests", "outputs", "give_completeness", output_filename
        )
        modules_table_file = os.path.join(
            repo_root, "tests", "fixtures", "give_completeness", "modules_table.tsv"
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        assert errors == []

    def test_combined_contigs_dataset(self):
        repo_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        output_file = os.path.join(
            repo_root, "tests", "fixtures", "utils", "sanity_check", "test_combined_contigs_broken.tsv"
        )
        modules_table_file = os.path.join(
            repo_root, "tests", "fixtures", "give_completeness", "modules_table.tsv"
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        errors_line = ';'.join(errors)
        assert "matching_ko K01810 not found in definition of module M00909" in errors_line
        assert "matching_ko K25026 not found in definition of module M00909" in errors_line
        assert "matching_ko K01810 not found in definition of module M00892" in errors_line
        assert "matching_ko K00844 not found in definition of module M00892" in errors_line
        assert "missing_ko K00844 not found in definition of module M00892" in errors_line
        assert "matching_ko K01810 not found in definition of module M00892" in errors_line
        assert "matching_ko K01810 not found in definition of module M00909" in errors_line
        assert "missing_ko K00845 not found in definition of module M00909" in errors_line
        assert len(errors) == 8

    def test_combined_pathways_dataset(self):
        repo_root = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        output_file = os.path.join(
            repo_root, "tests", "fixtures", "utils", "sanity_check", "test_combined_pathways_broken.tsv"
        )
        modules_table_file = os.path.join(
            repo_root, "tests", "fixtures", "give_completeness", "modules_table.tsv"
        )
        errors = sanity_check_give_completeness_output(output_file, modules_table_file)
        errors_line = ';'.join(errors)
        assert "matching_ko K01810 not found in definition of module M00909" in errors_line
        assert "matching_ko K25026 not found in definition of module M00909" in errors_line
        assert "matching_ko K00844 not found in definition of module M00892" in errors_line
        assert "matching_ko K01810 not found in definition of module M00892" in errors_line
        assert len(errors) == 4

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
