#!/bin/env python3

import unittest
import os
from hashlib import md5
from kegg_pathways_completeness.bin.fetch_modules_data import ModulesDataFetchTool


def calculate_md5(file_path):
    hash_md5 = md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


class TestFetchModulesData(unittest.TestCase):
    def test_fetch_module_info(self):
        output_dir = "output"
        module_fetch_tool = ModulesDataFetchTool(output_dir)
        module_fetch_tool.fetch_module_info("M00001")
        expected_files = {
            "all_pathways.txt": "8f34cd2ee36c25102a7ae621d2e014a1",
            "all_pathways_class.txt": "2e98048a59f8e72bc9feb1f5590f84a8",
            "all_pathways_names.txt": "b7813827ac438c2f255824c86ced4b81",
        }
        for file_name, expected_md5 in expected_files.items():
            file_path = os.path.join(output_dir, file_name)
            self.assertTrue(
                os.path.exists(file_path), f"File {file_name} does not exist"
            )
            actual_md5 = calculate_md5(file_path)
            self.assertEqual(
                actual_md5, expected_md5, f"MD5 mismatch for file {file_name}"
            )


if __name__ == "__main__":
    unittest.main()
