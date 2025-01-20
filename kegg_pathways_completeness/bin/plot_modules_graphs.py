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
import sys
import pickle
import networkx as nx
import logging
import os
import graphviz
import pydot
import csv

from .utils import parse_modules_list_input, parse_graphs_input, __version__


logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Script generates plots for each contig")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
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
    parser.add_argument("--use-pydot", dest="use_pydot", help="Use pydot instead of graphviz", required=False,
                        action='store_true')
    return parser.parse_args(argv)


class PlotModuleCompletenessGraph():
    def __init__(
            self,
            modules_completeness: dict,
            graphs: nx.MultiDiGraph,
            modules_definitions: dict,
            outdir: str,
            modules_list: list = [],
            use_pydot: bool = False
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
        self.modules_list = modules_list
        # flags
        self.use_pydot = use_pydot
        # output directories
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.outdir_dot = os.path.join(self.outdir, 'dot')
        if not os.path.exists(self.outdir_dot):
            os.mkdir(self.outdir_dot)
        self.outdir_png = os.path.join(self.outdir, 'png')
        if not os.path.exists(self.outdir_png):
            os.mkdir(self.outdir_png)

    def format_weight(self, weight):
        if weight in (0, 1):
            return str(weight)
        return f"1/{int(1 / weight)}"

    def create_graph_dot(self, name, presented, graph, pathways_schema):
        # Create a pydot graph
        dot = pydot.Dot(name=name, graph_type='digraph', comment=pathways_schema)
        edges = graph[0].edges

        for edge, count in zip(edges, range(len(edges))):
            from_node = edge[0]
            to_node = edge[1]
            number = edge[2]
            # Add nodes to the graph
            dot.add_node(pydot.Node(str(from_node)))
            dot.add_node(pydot.Node(str(to_node)))
            # Extract edge attributes
            label = edges._adjdict[from_node][to_node][number]['label']
            weight = edges._adjdict[from_node][to_node][number]['weight']
            # Set edge color
            color = 'red' if label in presented else 'black'
            # Add edge to the graph
            dot.add_edge(
                pydot.Edge(
                    str(from_node),
                    str(to_node),
                    label=f"{label} \n [{self.format_weight(weight)}]",
                    color=color
                )
            )
        return dot

    def create_graph(self, name, presented, graph, pathways_schema):
        dot = graphviz.Digraph(name, comment=pathways_schema)
        edges = graph[0].edges
        for edge, count in zip(edges, range(len(edges))):
            from_node = edge[0]
            to_node = edge[1]
            number = edge[2]
            dot.node(str(from_node))
            dot.node(str(to_node))

            label = edges._adjdict[from_node][to_node][number]['label']
            weight = edges._adjdict[from_node][to_node][number]['weight']

            color = 'red' if label in presented else 'black'
            dot.edge(str(from_node), str(to_node), label=label + ' \n [' + self.format_weight(weight) + ']', color=color)
        return dot

    def generate_graph_using_pydot(self, name, presented_ko, graph, pathways_schema):
        logging.info('Using pydot')
        dot = self.create_graph_dot(name, presented=presented_ko, graph=graph,
                                    pathways_schema=pathways_schema)
        # create .dot file
        with open(os.path.join(self.outdir_dot, f"{name}.dot"), "w") as f:
            f.write(dot.to_string())
        # create .png file
        dot.write_png(os.path.join(self.outdir_png, f'{name}.png'))

    def generate_graph_using_graphviz(self, name, presented_ko, graph, pathways_schema):
        logging.info('Using graphviz')
        dot = self.create_graph(name, presented=presented_ko, graph=graph,
                                pathways_schema=pathways_schema)
        dot.render(directory=self.outdir, filename=name, format='png')

    def generate_plot_for_completeness(self):
        logging.info('Using completeness file')
        for name in self.modules_completeness:
            if len(self.modules_list):
                if name not in self.modules_list:
                    logging.debug(f'Skipping {name} because it is not in specified modules list')
                    continue
            graph = self.graphs[name]
            logging.info(f'Plotting {name}')
            presented_ko = self.modules_completeness[name].split(',')
            if self.use_pydot:
                self.generate_graph_using_pydot(name=name, presented_ko=presented_ko, graph=graph,
                                                pathways_schema=self.modules_definitions[name])
            else:
                self.generate_graph_using_graphviz(name=name, presented_ko=presented_ko, graph=graph,
                                                   pathways_schema=self.modules_definitions[name])

    def generate_plot_without_completeness(self):
        logging.info("Plotting modules from specified list without completeness information")
        for name in self.modules_list:
            graph = self.graphs[name]
            logging.info(f'Plotting {name}')
            if self.use_pydot:
                self.generate_graph_using_pydot(name=name, presented_ko=[], graph=graph,
                                                pathways_schema=self.modules_definitions[name])
            else:
                self.generate_graph_using_graphviz(name=name, presented_ko=[], graph=graph,
                                                   pathways_schema=self.modules_definitions[name])

    def generate_plot(self):
        if self.modules_completeness:
            self.generate_plot_for_completeness()
        else:
            self.generate_plot_without_completeness()


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


def parse_input_modules(input_modules_list=None, input_modules_file=None, list_separator='\n'):
    if input_modules_list:
        logging.info("Using specified list of modules")
        return input_modules_list
    elif input_modules_file:
        logging.info(f"Using modules from {input_modules_file}")
        modules_list = []
        with open(input_modules_file, 'r') as f:
            reader = csv.reader(f, delimiter=list_separator)
            for row in reader:
                modules_list.extend(row)
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
        modules_list=parse_input_modules(args.input_modules_list, args.input_modules_file, args.list_separator),
        use_pydot=args.use_pydot
    )
    plot_completeness_generator.generate_plot()


if __name__ == "__main__":
    main()

