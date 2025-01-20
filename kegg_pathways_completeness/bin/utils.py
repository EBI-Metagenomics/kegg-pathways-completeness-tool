#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import logging
import networkx as nx
import os
from pathlib import Path
import pickle
import sys
import tomli

__version__ = "1.2.1"

def setup_logging(verbose):
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s %(levelname)s - %(message)s'
    )


def intersection(lst1, lst2):
    """
    Intersection of two sets
    :param lst1: first input list
    :param lst2: second input list
    :return: intersection in list format
    """
    return list(set(lst1) & set(lst2))


def parse_modules_list_input(filepath):
    """
    Function parses input file with modules definitions to generate dictionary in format:
    [module]:KOs schema
    :param filepath: all_pathways.txt
    :return: dictionary [module]:KOs schema
    """
    pathways_schema = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as pathways_file:
            for line in pathways_file:
                fields = line.strip().split(':')
                pathways_schema[fields[0]] = fields[1]
        return pathways_schema
    else:
        logging.error(f'No file {filepath} found')


def parse_graphs_input(filename):
    """
    Function loads graphs of modules in networkx format pre-saved into pkl format.
    :param filename: graphs.pkl
    :return: Graph
    """
    if os.path.exists(filename):
        with open(filename, 'rb') as file_graph:
            graphs = pickle.load(file_graph)
        return graphs
    else:
        logging.error(f'No graphs {filename} file found')