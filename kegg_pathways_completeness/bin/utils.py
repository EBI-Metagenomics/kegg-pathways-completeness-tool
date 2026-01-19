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
