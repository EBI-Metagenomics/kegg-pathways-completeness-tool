#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2026 EMBL - European Bioinformatics Institute
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


import logging
import os
import pickle
import re
from importlib.metadata import PackageNotFoundError, version


def get_version():
    """Get package version from installed metadata"""
    try:
        return version("kegg-pathways-completeness")
    except PackageNotFoundError:
        return "unknown"


def setup_logging(verbose):
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
    )


def intersection(lst1, lst2):
    """
    Intersection of two sets
    :param lst1: first input list
    :param lst2: second input list
    :return: intersection in list format
    """
    return list(set(lst1) & set(lst2))


def parse_graphs_input(filename):
    """
    Function loads graphs of modules in networkx format pre-saved into pkl format.
    :param filename: graphs.pkl
    :return: Graph
    """
    if os.path.exists(filename):
        with open(filename, "rb") as file_graph:
            graphs = pickle.load(file_graph)
        return graphs
    else:
        logging.error(f"No graphs {filename} file found")


def parse_module_kos(definition):
    """
    Extract the set of KO accessions referenced in a module definition string.
    :param definition: module definition, ex. "(K00844,K12407) K01803"
    :return: set of KO accessions, ex. {"K00844", "K12407", "K01803"}
    """
    return set(re.findall(r"K\d{5}", definition))


def load_modules_kos(modules_table_file):
    """
    Load KO accessions referenced in each module's definition from modules_table.tsv.
    :param modules_table_file: modules_table.tsv (columns: module, definition, name, class)
    :return: dict {module_accession: set of KOs in that module's definition}
    """
    modules_kos = {}
    with open(modules_table_file, "r") as f:
        header = f.readline().strip().split("\t")
        for col in ("module", "definition"):
            if col not in header:
                raise ValueError(f"TSV file must have '{col}' column")
        module_idx = header.index("module")
        definition_idx = header.index("definition")

        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > max(module_idx, definition_idx):
                module = fields[module_idx]
                modules_kos[module] = parse_module_kos(fields[definition_idx])

    return modules_kos


def sanity_check_give_completeness_output(output_file, modules_table_file):
    """
    Sanity check a give_completeness output table (pathways or per-contig summary):
    for each row, both matching_ko and missing_ko must exist in that module's
    definition (both columns list KOs that are part of the module's structure -
    matching_ko are present in the input, missing_ko are required but absent).
    :param output_file: give_completeness output TSV
    (columns include module_accession, matching_ko, missing_ko)
    :param modules_table_file: modules_table.tsv used to run give_completeness
    :return: list of error messages (empty list means the output is consistent)
    """
    logger = logging.getLogger(__name__)
    modules_kos = load_modules_kos(modules_table_file)

    errors = []
    with open(output_file, "r") as f:
        header = f.readline().strip().split("\t")
        for col in ("module_accession", "matching_ko", "missing_ko"):
            if col not in header:
                raise ValueError(f"TSV file must have '{col}' column")
        module_idx = header.index("module_accession")
        matching_idx = header.index("matching_ko")
        missing_idx = header.index("missing_ko")

        for line_num, line in enumerate(f, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            module = fields[module_idx]

            if module not in modules_kos:
                errors.append(
                    f"Line {line_num}: module {module} not found in {modules_table_file}"
                )
                continue
            definition_kos = modules_kos[module]

            matching_ko = (
                parse_module_kos(fields[matching_idx])
                if len(fields) > matching_idx
                else set()
            )
            missing_ko = (
                parse_module_kos(fields[missing_idx])
                if len(fields) > missing_idx
                else set()
            )

            for ko in sorted(matching_ko - definition_kos):
                errors.append(
                    f"Line {line_num}: matching_ko {ko} not found in "
                    f"definition of module {module}"
                )
            for ko in sorted(missing_ko - definition_kos):
                errors.append(
                    f"Line {line_num}: missing_ko {ko} not found in "
                    f"definition of module {module}"
                )

    if errors:
        for error in errors:
            logger.error(error)
    else:
        logger.info(f"Sanity check passed for {output_file}")

    return errors
