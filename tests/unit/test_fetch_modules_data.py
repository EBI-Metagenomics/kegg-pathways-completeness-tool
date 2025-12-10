#!/usr/bin/env python3

import os
import pytest
import tempfile
import shutil
from unittest.mock import Mock, patch, MagicMock
from kegg_pathways_completeness.bin.fetch_modules_data import (
    ModulesDataFetchTool,
    get_version,
)


# Mock KEGG API responses
MOCK_MODULE_LIST_RESPONSE = """md:M00001\tGlycolysis
md:M00002\tCitrate cycle
md:M00003\tGluconneogenesis"""

MOCK_MODULE_M00001_RESPONSE = """ENTRY       M00001            Pathway   Module
NAME        Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
DEFINITION  K00844 K12407 K00845 K01810
CLASS       Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism
"""

MOCK_MODULE_M00002_RESPONSE = """ENTRY       M00002            Pathway   Module
NAME        Citrate cycle (TCA cycle), aerobic respiration
DEFINITION  K01647 K01681 K01682 K01902
CLASS       Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism
"""

MOCK_MODULE_M00003_RESPONSE = """ENTRY       M00003            Pathway   Module
NAME        Gluconeogenesis, oxaloacetate => fructose-6P
DEFINITION  K01596 K01610 K01834
CLASS       Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism
"""


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs"""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def mock_session():
    """Create a mock session with predefined responses"""
    with patch(
        "kegg_pathways_completeness.bin.fetch_modules_data.requests.Session"
    ) as mock:
        session = MagicMock()
        mock.return_value = session

        def get_side_effect(url, timeout=30):
            response = Mock()
            response.raise_for_status = Mock()

            if "list/module" in url:
                response.text = MOCK_MODULE_LIST_RESPONSE
            elif "M00001" in url:
                response.text = MOCK_MODULE_M00001_RESPONSE
            elif "M00002" in url:
                response.text = MOCK_MODULE_M00002_RESPONSE
            elif "M00003" in url:
                response.text = MOCK_MODULE_M00003_RESPONSE
            else:
                response.raise_for_status.side_effect = Exception("Unknown URL")

            return response

        session.get.side_effect = get_side_effect
        yield session


class TestModulesDataFetchTool:
    """Test suite for ModulesDataFetchTool"""

    def test_init(self, temp_dir):
        """Test tool initialization"""
        tool = ModulesDataFetchTool(
            outdir=temp_dir, max_workers=5, max_retries=3, delay_between_requests=0.1
        )
        assert tool.outdir == temp_dir
        assert tool.max_workers == 5
        assert tool.max_retries == 3
        assert tool.delay_between_requests == 0.1
        assert tool.modules_data == []
        assert tool.failed_modules == []
        assert os.path.exists(temp_dir)

    def test_fetch_and_save_kegg_modules_list(self, temp_dir, mock_session):
        """Test fetching module list from API"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        modules = tool.fetch_and_save_kegg_modules_list()

        assert len(modules) == 3
        assert "md:M00001" in modules
        assert "md:M00002" in modules
        assert "md:M00003" in modules

        # Check file was created
        list_file = os.path.join(temp_dir, "modules_list.txt")
        assert os.path.exists(list_file)

        with open(list_file, "r") as f:
            content = f.read()
            assert "md:M00001" in content
            assert "md:M00002" in content

    @patch("kegg_pathways_completeness.bin.fetch_modules_data.time.sleep")
    def test_fetch_module_info(self, mock_sleep, temp_dir, mock_session):
        """Test fetching single module information"""
        tool = ModulesDataFetchTool(outdir=temp_dir, delay_between_requests=0.1)
        result = tool.fetch_module_info("M00001")

        assert result is not None
        assert result["module"] == "M00001"
        assert "Glycolysis" in result["name"]
        assert "K00844" in result["definition"]
        assert "Carbohydrate metabolism" in result["class"]

        # Verify sleep was called for rate limiting
        mock_sleep.assert_called_with(0.1)

    @patch("kegg_pathways_completeness.bin.fetch_modules_data.time.sleep")
    def test_fetch_all_modules_parallel(self, mock_sleep, temp_dir, mock_session):
        """Test parallel fetching of multiple modules"""
        tool = ModulesDataFetchTool(
            outdir=temp_dir, max_workers=2, delay_between_requests=0.0
        )
        modules = ["M00001", "M00002", "M00003"]

        tool.fetch_all_modules_parallel(modules)

        assert len(tool.modules_data) == 3
        assert all(
            m["module"] in ["M00001", "M00002", "M00003"] for m in tool.modules_data
        )

    def test_write_modules_table(self, temp_dir):
        """Test writing modules table to TSV file"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        tool.modules_data = [
            {
                "module": "M00001",
                "definition": "K00844 K12407",
                "name": "Glycolysis",
                "class": "Carbohydrate metabolism",
            },
            {
                "module": "M00002",
                "definition": "K01647 K01681",
                "name": "TCA cycle",
                "class": "Carbohydrate metabolism",
            },
        ]

        tool.write_modules_table()

        table_file = os.path.join(temp_dir, "modules_table.tsv")
        assert os.path.exists(table_file)

        with open(table_file, "r") as f:
            lines = f.readlines()
            assert len(lines) == 3  # header + 2 data rows
            assert lines[0].strip() == "module\tdefinition\tname\tclass"
            assert "M00001" in lines[1]
            assert "M00002" in lines[2]

    def test_parse_old_file(self, temp_dir):
        """Test parsing old format files"""
        old_file = os.path.join(temp_dir, "old_definitions.txt")
        with open(old_file, "w") as f:
            f.write("M00001:K00844 K12407\n")
            f.write("M00002:K01647 K01681\n")

        result = ModulesDataFetchTool.parse_old_file(old_file)

        assert len(result) == 2
        assert result["M00001"] == "K00844 K12407"
        assert result["M00002"] == "K01647 K01681"

    def test_detect_changes_new_module(self, temp_dir):
        """Test detecting new modules"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        tool.modules_data = [
            {
                "module": "M00001",
                "definition": "K00844 K12407",
                "name": "Glycolysis",
                "class": "Carbohydrate metabolism",
            },
            {
                "module": "M00003",  # New module
                "definition": "K01596 K01610",
                "name": "Gluconeogenesis",
                "class": "Carbohydrate metabolism",
            },
        ]

        old_definitions = {
            "M00001": "K00844 K12407"
            # M00003 is missing (new module)
        }

        changed = tool.detect_changes(old_definitions=old_definitions)

        assert len(changed) == 1
        assert changed[0]["module"] == "M00003"

    def test_detect_changes_modified_definition(self, temp_dir):
        """Test detecting modified definitions"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        tool.modules_data = [
            {
                "module": "M00001",
                "definition": "K00844 K12407 K00845",  # Modified
                "name": "Glycolysis",
                "class": "Carbohydrate metabolism",
            }
        ]

        old_definitions = {"M00001": "K00844 K12407"}  # Old definition

        changed = tool.detect_changes(old_definitions=old_definitions)

        assert len(changed) == 1
        assert changed[0]["module"] == "M00001"

    def test_detect_changes_modified_name(self, temp_dir):
        """Test detecting modified names"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        tool.modules_data = [
            {
                "module": "M00001",
                "definition": "K00844 K12407",
                "name": "Glycolysis (updated)",  # Modified name
                "class": "Carbohydrate metabolism",
            }
        ]

        old_names = {"M00001": "Glycolysis"}

        changed = tool.detect_changes(old_names=old_names)

        assert len(changed) == 1
        assert changed[0]["module"] == "M00001"

    def test_detect_no_changes(self, temp_dir):
        """Test when there are no changes"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        tool.modules_data = [
            {
                "module": "M00001",
                "definition": "K00844 K12407",
                "name": "Glycolysis",
                "class": "Carbohydrate metabolism",
            }
        ]

        old_definitions = {"M00001": "K00844 K12407"}
        old_names = {"M00001": "Glycolysis"}
        old_classes = {"M00001": "Carbohydrate metabolism"}

        changed = tool.detect_changes(
            old_definitions=old_definitions,
            old_names=old_names,
            old_classes=old_classes,
        )

        assert len(changed) == 0

    def test_write_changed_table(self, temp_dir):
        """Test writing changed modules table"""
        tool = ModulesDataFetchTool(outdir=temp_dir)
        changed_modules = [
            {
                "module": "M00001",
                "definition": "K00844 K12407 K00845",
                "name": "Glycolysis (updated)",
                "class": "Carbohydrate metabolism",
            }
        ]

        tool.write_changed_table(changed_modules)

        changed_file = os.path.join(temp_dir, "changed.tsv")
        assert os.path.exists(changed_file)

        with open(changed_file, "r") as f:
            lines = f.readlines()
            assert len(lines) == 2  # header + 1 data row
            assert "M00001" in lines[1]
            assert "updated" in lines[1]

    def test_get_version(self):
        """Test version retrieval"""
        version = get_version()
        assert version is not None
        # Version should be either a valid version string or 'unknown'
        assert isinstance(version, str)


class TestRetryMechanism:
    """Test suite for retry and rate limiting features"""

    @patch("kegg_pathways_completeness.bin.fetch_modules_data.time.sleep")
    def test_rate_limiting_delay(self, mock_sleep, temp_dir, mock_session):
        """Test that delay is applied between requests"""
        tool = ModulesDataFetchTool(outdir=temp_dir, delay_between_requests=0.5)
        tool.fetch_module_info("M00001")

        mock_sleep.assert_called_with(0.5)

    def test_failed_module_tracking(self, temp_dir):
        """Test tracking of failed module fetches"""
        tool = ModulesDataFetchTool(outdir=temp_dir)

        with patch.object(tool.session, "get") as mock_get:
            mock_get.side_effect = Exception("Network error")
            result = tool.fetch_module_info("M00999")

            assert result is None
            assert "M00999" in tool.failed_modules


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
