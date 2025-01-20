#!/usr/bin/env python3
# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os
import subprocess
import sys
import re
from .utils import __version__


LIST_MODULES = 'list_pathways.txt'
LIST_CLASSES = "all_pathways_class.txt"
LIST_NAMES = "all_pathways_names.txt"
LIST_PATHS = "all_pathways.txt"
LIST_SEPARATED_IN_DEFINITION = "definition_separated.txt"

KEGG_API_LIST_MODULES = "http://rest.kegg.jp/list/module"
KEGG_API_GET = "http://rest.kegg.jp/get"


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Script fetches KEGG API for list of modules with NAME, DEFINITION and CLASS.")
    parser.add_argument("-o", "--output-dir", dest="output", help="Output directory", default="pathways",
                        required=False)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args(argv)


class ModulesDataFetchTool():
    def __init__(
            self,
            outdir: str,
    ):
        """
        This script uses subprocess to get list of KEGG modules via KEGG API.
        By default it generates tree files:
        all_pathways.txt - list of all modules in format module:kos schema (MXXXXX:KYYYYY+KZZZZZ)
        all_pathways_class.txt - list of CLASS definitions of modules in format module:class
        all_pathways_names.txt - list of NAME definitions of modules in format module:name
        It is useful to run whole fetch and update because KEGG changes modules structure for already existing modules.
        :param outdir: path to output directory
        :return:
        """
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    def fetch_and_save_kegg_modules_list(self):
        # Ensure the output directory exists
        os.makedirs(self.outdir, exist_ok=True)
        output_file = os.path.join(self.outdir, LIST_MODULES)

        try:
            print(f'Fetching modules from {KEGG_API_LIST_MODULES}')
            # Fetch the data from the KEGG API using wget
            wget_command = f'wget -qO- {KEGG_API_LIST_MODULES}'
            wget_result = subprocess.run(wget_command, shell=True, capture_output=True, text=True)
            wget_output = wget_result.stdout

            # Process the output to extract the module codes
            module_codes = [line.split('\t')[0] for line in wget_output.splitlines()]
            # Write the results to the output file
            with open(output_file, 'w') as f:
                for code in module_codes:
                    f.write(code + '\n')

            print(f'Successfully saved KEGG module codes to {output_file}')
            return module_codes
        except subprocess.CalledProcessError as e:
            print(f'An error occurred: {e}')
            exit(1)

    def fetch_module_info(self, module):
        try:
            line_separator_in_kegg_names = ' '  # space=AND, comma=OR
            print(f'Fetching modules from {KEGG_API_GET}/{module}')
            # Fetch the data from the KEGG API using wget
            wget_command = f'wget -qO- {KEGG_API_GET}/{module}'
            wget_result = subprocess.run(wget_command, shell=True, capture_output=True, text=True)
            wget_output = wget_result.stdout
            lines = wget_output.split('\n')
            definition_lines = []
            capture_definition = False
            for line in lines:
                if 'NAME' in line:
                    name = re.search(r'^NAME\s+(.*)', line).group(1)
                elif 'CLASS' in line:
                    class_module = re.search(r'^CLASS\s+(.*)', line).group(1)
                elif 'DEFINITION' in line:
                    # can be in multiple lines
                    definition_lines.append(re.search(r'^DEFINITION\s+(.*)', line).group(1))
                    capture_definition = True
                elif capture_definition and line.startswith(" "):
                    definition_lines.append(line.strip())
                else:
                    capture_definition = False
            if len(definition_lines) > 1:
                for i in range(len(definition_lines)):
                    if len(definition_lines[i]) > 6:
                        # not a single KOXXXX
                        definition_lines[i] = f"({definition_lines[i]})"
            def_module = line_separator_in_kegg_names.join(definition_lines)
            with open(os.path.join(self.outdir, LIST_PATHS), 'a') as file_out:
                file_out.write(f'{module}:{def_module}' + '\n')
            with open(os.path.join(self.outdir, LIST_CLASSES), 'a') as file_out:
                file_out.write(f'{module}:{class_module}' + '\n')
            with open(os.path.join(self.outdir, LIST_NAMES), 'a') as file_out:
                file_out.write(f'{module}:{name}' + '\n')
            if len(definition_lines) > 1:
                with open(os.path.join(self.outdir, LIST_SEPARATED_IN_DEFINITION), 'a') as file_out:
                    file_out.write(f'{module}:{def_module}' + '\n')

        except subprocess.CalledProcessError as e:
            print(f'An error occurred: {e}')
            exit(1)


def main():
    args = parse_args(sys.argv[1:])
    module_fetch_tool = ModulesDataFetchTool(
        outdir=args.output
    )
    # get a list of all modules available in KEGG API
    modules_all = module_fetch_tool.fetch_and_save_kegg_modules_list()
    # get info for each module
    for module in modules_all:
        module_fetch_tool.fetch_module_info(module)


if __name__ == "__main__":
    main()
