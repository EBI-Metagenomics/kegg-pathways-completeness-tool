#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
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
import sys
import pickle
import networkx as nx
import logging
import os
import graphviz

from .utils import parse_modules_list_input, parse_graphs_input

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Script generates plots for each contig")
    parser.add_argument("-i", "--input-completeness", dest="input_completeness",
                        help="Output table from give_completeness.py", required=False)
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-m", "--modules", dest="input_modules_list", nargs='+',
                       help="Space separated list of modules accessions")
    group.add_argument("-l", "--modules-file", dest="input_modules_file",
                       help="File containing modules accessions")
    group.add_argument("-s", "--file-separator", dest="list_separator",
                       help="Modules separator in file", default='\n')
    parser.add_argument("-g", "--graphs", dest="graphs", help="graphs in pickle format", required=False,
                        default="pathways_data/graphs.pkl")
    parser.add_argument("-d", "--definitions", dest="pathways", help="Pathways of kos", required=False,
                        default="pathways_data/all_pathways.txt")
    parser.add_argument("-o", "--outdir", dest="outdir", help="Path to output directory", required=False,
                        default="pathways_plots")
    return parser.parse_args(argv)


class PlotModuleCompletenessGraph():
    def __init__(
            self,
            modules_completeness: dict,
            graphs: nx.MultiDiGraph,
            modules_definitions: dict,
            outdir: str,
            modules_list: list = [],
    ):
        """
        Class generates a plots with colored edges depending on presence/absence.
        If edge is presented - it is colored in red, otherwise - black
        :param modules_completeness: parsed output of give_completeness.py
        :param graphs: graphs in nx format
        :param modules_definitions: pre-parsed all_pathways.txt
        :param outdir: name of output directory, default=pathways_plots
        :param modules_list: list of modules of interest
        """
        self.modules_completeness = modules_completeness
        self.graphs = graphs
        self.modules_definitions = modules_definitions
        self.outdir = outdir
        self.modules_list = modules_list

    def create_graph(self, name, presented, graph, pathways_schema):
        dot = graphviz.Digraph(name, comment=pathways_schema)
        edges = graph[0].edges
        max_weight = 0
        for edge, count in zip(edges, range(len(edges))):
            from_node = edge[0]
            to_node = edge[1]
            number = edge[2]
            dot.node(str(from_node))
            dot.node(str(to_node))

            label = edges._adjdict[from_node][to_node][number]['label']
            weight = edges._adjdict[from_node][to_node][number]['weight']
            if weight > 0:
                if 1 / weight > max_weight:
                    max_weight = int(1 / weight)
            if weight == 1 or weight == 0:
                weight_str = str(weight)
            else:
                weight_str = '1/' + str(int(1 / weight))
            color = 'red' if label in presented else 'black'
            dot.edge(str(from_node), str(to_node), label=label + ' \n [' + weight_str + ']', color=color)
        return dot

    def generate_plot_for_completeness(self):
        if self.modules_completeness:
            logging.info('Using completeness file')
            for name in self.modules_completeness:
                if len(self.modules_list):
                    if name not in self.modules_list:
                        logging.debug(f'Skipping {name} because it is not in specified modules list')
                        continue
                graph = self.graphs[name]
                logging.info(f'Plotting {name}')
                presented_ko = self.modules_completeness[name].split(',')
                dot = self.create_graph(name, presented=presented_ko, graph=graph,
                                        pathways_schema=self.modules_definitions[name])
                dot.render(directory=self.outdir, filename=name, format='png')
        else:
            logging.info("Plotting modules from specified list without completeness information")
            for name in self.modules_list:
                graph = self.graphs[name]
                logging.info(f'Plotting {name}')
                dot = self.create_graph(name, presented=[], graph=graph,
                                        pathways_schema=self.modules_definitions[name])
                dot.render(directory=self.outdir, filename=name, format='png')


def parse_completeness_input(filepath):
    """
    Function parses input file with completeness to generate dictionary in format:
    [module]:matching_ko
    """
    pathways = {}
    if filepath:
        with open(filepath, 'r') as file_in:
            for line in file_in:
                line = line.strip().split('\t')
                if line[0] == 'contig':
                    column_pathway = 1
                    column_kos = 5
                elif line[0] == 'module_accession':
                    column_pathway = 0
                    column_kos = 4
                else:
                    pathways[line[column_pathway]] = line[column_kos]
    return pathways


def parse_input_modules(input_modules_list=None, input_modules_file=None):
    if input_modules_list:
        logging.info("Using specified list of modules")
        return input_modules_list
    elif input_modules_file:
        logging.info(f"Using modules from {input_modules_file}")
        modules_list = []
        with open(input_modules_file, 'r') as file_in:
            for line in file_in:
                modules_list.append(line.strip())
                # TODO add separator
        return modules_list
    else:
        logging.info('No modules specified in input')
        return []


def main():
    args = parse_args(sys.argv[1:])
    if not args.input_completeness and not args.input_modules_list and not args.input_modules_file:
        logging.error("None of input files was presented. Please, add input file with -i or -m or -l")
        exit(1)
    plot_completeness_generator = PlotModuleCompletenessGraph(
        modules_completeness=parse_completeness_input(args.input_completeness),
        graphs=parse_graphs_input(args.graphs),
        modules_definitions=parse_modules_list_input(args.pathways),
        outdir=args.outdir,
        modules_list=parse_input_modules(args.input_modules_list, args.input_modules_file),
    )
    plot_completeness_generator.generate_plot_for_completeness()


if __name__ == "__main__":
    main()
