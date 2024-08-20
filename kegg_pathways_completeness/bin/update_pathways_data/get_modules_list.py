#!/usr/bin/env python3
"""
This script uses subprocess to get list of KEGG modules via KEGG API.

By default it generates tree files:
all_pathways.txt - list of all modules in format module:kos schema (MXXXXX:KYYYYY+KZZZZZ)
all_pathways_class.txt - list of CLASS definitions of modules in format module:class
all_pathways_names.txt - list of NAME definitions of modules in format module:name
It is useful to run whole fetch and update because KEGG changes modules structure for already existing modules.

It is also possible to check only what new modules have appeared in comparison with previous version in this repository.
Script fetches a list of modules and compares it with pathways_data/all_pathways.txt file.
It produces file new_modules.txt with list of missing modules.
"""

import argparse
import os
import subprocess
import sys
import re


LIST_CLASSES = "all_pathways_class.txt"
LIST_NAMES = "all_pathways_names.txt"
LIST_PATHS = "all_pathways.txt"

def parse_args():
    parser = argparse.ArgumentParser(
        description="Script fetches KEGG API for list of modules with NAME, DEFINITION and CLASS.")
    parser.add_argument("--check-new-modules", dest="check_new_modules", help="Check only what new modules appeared",
                        action='store_true')
    parser.add_argument("-o", "--output-dir", dest="output", help="Output directory", default="pathways")
    return parser


def fetch_and_save_kegg_modules_list(output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'list_pathways.txt')

    try:
        print(f'Fetching modules from http://rest.kegg.jp/list/module')
        # Fetch the data from the KEGG API using wget
        wget_command = f'wget -qO- http://rest.kegg.jp/list/module'
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


def compare_with_existing(new_modules):
    existing = []
    with open('kegg_pathways_completeness/pathways_data/all_pathways.txt', 'r') as file_in:
        for line in file_in:
            existing.append(line.strip().split(':')[0])
    result = list(set(new_modules).difference(set(existing)))
    if result:
        print(f'Found {len(result)} new modules')
        with open('new_modules.txt', 'w') as file_out:
            file_out.write('\n'.join(result))
        return result
    else:
        print('No new modules')
        return False


def fetch_module_info(module, output_dir):
    try:
        print(f'Fetching modules from http://rest.kegg.jp/get/{module}')
        # Fetch the data from the KEGG API using wget
        wget_command = f'wget -qO- http://rest.kegg.jp/get/{module}'
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
        def_module = ' '.join(definition_lines)
        with open(os.path.join(output_dir, LIST_PATHS), 'a') as file_out:
            file_out.write(f'{module}:{def_module}' + '\n')
        with open(os.path.join(output_dir, LIST_CLASSES), 'a') as file_out:
            file_out.write(f'{module}:{class_module}' + '\n')
        with open(os.path.join(output_dir, LIST_NAMES), 'a') as file_out:
            file_out.write(f'{module}:{name}' + '\n')

    except subprocess.CalledProcessError as e:
        print(f'An error occurred: {e}')
        exit(1)


if __name__ == '__main__':
    parser = parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        # get a list of all modules available in KEGG API
        modules_all = fetch_and_save_kegg_modules_list(args.output)
        if args.check_new_modules:
            # find missing modules
            new_modules = compare_with_existing(modules_all)
        else:
            # regenerate for all modules
            new_modules = modules_all
        if new_modules:
            # new modules found
            for module in new_modules:
                fetch_module_info(module, args.output)

