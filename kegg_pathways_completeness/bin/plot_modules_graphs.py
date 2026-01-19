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


import csv
import logging
import os
from importlib.resources import files

import click
import graphviz
import networkx as nx
import pydot

from .utils import get_version, parse_graphs_input

logging.basicConfig(encoding="utf-8", level=logging.DEBUG)


class PlotModuleCompletenessGraph:
    def __init__(
        self,
        modules_completeness: dict,
        graphs: nx.MultiDiGraph,
        modules_definitions: dict,
        outdir: str,
        modules_list: list = [],
        use_pydot: bool = False,
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
        self.outdir_dot = os.path.join(self.outdir, "dot")
        if not os.path.exists(self.outdir_dot):
            os.mkdir(self.outdir_dot)
        self.outdir_png = os.path.join(self.outdir, "png")
        if not os.path.exists(self.outdir_png):
            os.mkdir(self.outdir_png)

    def format_weight(self, weight):
        if weight in (0, 1):
            return str(weight)
        return f"1/{int(1 / weight)}"

    def create_graph_dot(self, name, presented, graph, pathways_schema):
        # Create a pydot graph
        dot = pydot.Dot(name=name, graph_type="digraph", comment=pathways_schema)
        edges = graph[0].edges

        for edge, count in zip(edges, range(len(edges))):
            from_node = edge[0]
            to_node = edge[1]
            number = edge[2]
            # Add nodes to the graph
            dot.add_node(pydot.Node(str(from_node)))
            dot.add_node(pydot.Node(str(to_node)))
            # Extract edge attributes
            label = edges._adjdict[from_node][to_node][number]["label"]
            weight = edges._adjdict[from_node][to_node][number]["weight"]
            # Set edge color
            color = "red" if label in presented else "black"
            # Add edge to the graph
            dot.add_edge(
                pydot.Edge(
                    str(from_node),
                    str(to_node),
                    label=f"{label} \n [{self.format_weight(weight)}]",
                    color=color,
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

            label = edges._adjdict[from_node][to_node][number]["label"]
            weight = edges._adjdict[from_node][to_node][number]["weight"]

            color = "red" if label in presented else "black"
            dot.edge(
                str(from_node),
                str(to_node),
                label=label + " \n [" + self.format_weight(weight) + "]",
                color=color,
            )
        return dot

    def generate_graph_using_pydot(self, name, presented_ko, graph, pathways_schema):
        logging.info("Using pydot")
        dot = self.create_graph_dot(
            name, presented=presented_ko, graph=graph, pathways_schema=pathways_schema
        )
        # create .dot file
        with open(os.path.join(self.outdir_dot, f"{name}.dot"), "w") as f:
            f.write(dot.to_string())
        # create .png file
        dot.write_png(os.path.join(self.outdir_png, f"{name}.png"))

    def generate_graph_using_graphviz(self, name, presented_ko, graph, pathways_schema):
        logging.info("Using graphviz")
        dot = self.create_graph(
            name, presented=presented_ko, graph=graph, pathways_schema=pathways_schema
        )
        dot.render(directory=self.outdir, filename=name, format="png")

    def generate_plot_for_completeness(self):
        logging.info("Using completeness file")
        for name in self.modules_completeness:
            if len(self.modules_list):
                if name not in self.modules_list:
                    logging.debug(
                        f"Skipping {name} because it is not in specified modules list"
                    )
                    continue
            graph = self.graphs[name]
            logging.info(f"Plotting {name}")
            presented_ko = self.modules_completeness[name].split(",")
            pathways_schema = (
                self.modules_definitions[name] if self.modules_definitions else name
            )
            if self.use_pydot:
                self.generate_graph_using_pydot(
                    name=name,
                    presented_ko=presented_ko,
                    graph=graph,
                    pathways_schema=pathways_schema,
                )
            else:
                self.generate_graph_using_graphviz(
                    name=name,
                    presented_ko=presented_ko,
                    graph=graph,
                    pathways_schema=pathways_schema,
                )

    def generate_plot_without_completeness(self):
        logging.info(
            "Plotting modules from specified list without completeness information"
        )
        for name in self.modules_list:
            graph = self.graphs[name]
            logging.info(f"Plotting {name}")
            pathways_schema = (
                self.modules_definitions[name] if self.modules_definitions else name
            )
            if self.use_pydot:
                self.generate_graph_using_pydot(
                    name=name,
                    presented_ko=[],
                    graph=graph,
                    pathways_schema=pathways_schema,
                )
            else:
                self.generate_graph_using_graphviz(
                    name=name,
                    presented_ko=[],
                    graph=graph,
                    pathways_schema=pathways_schema,
                )

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
        with open(filepath, "r") as file_in:
            for line in file_in:
                line = line.strip().split("\t")
                if line[0] == "contig":
                    column_pathway = 1
                    column_kos = 5
                elif line[0] == "module_accession":
                    column_pathway = 0
                    column_kos = 4
                else:
                    pathways[line[column_pathway]] = line[column_kos]
    return pathways


def parse_input_modules(
    input_modules_list=None, input_modules_file=None, list_separator="\n"
):
    if input_modules_list:
        logging.info("Using specified list of modules")
        return input_modules_list
    elif input_modules_file:
        logging.info(f"Using modules from {input_modules_file}")
        modules_list = []
        with open(input_modules_file, "r") as f:
            reader = csv.reader(f, delimiter=list_separator)
            for row in reader:
                modules_list.extend(row)
        return modules_list
    else:
        logging.info("No modules specified in input")
        return []


@click.command()
@click.option(
    "-i",
    "--input-completeness",
    type=click.Path(exists=True),
    help="Output table from give_completeness.py",
)
@click.option(
    "-m",
    "--modules",
    "input_modules_list",
    multiple=True,
    help="Module accessions (can be specified multiple times)",
)
@click.option(
    "-l",
    "--modules-file",
    type=click.Path(exists=True),
    help="File containing module accessions",
)
@click.option(
    "-s",
    "--file-separator",
    default="\n",
    help="Modules separator in file",
    show_default=True,
)
@click.option(
    "-g",
    "--graphs",
    type=click.Path(exists=True),
    help="Graphs in pickle format (default: uses packaged graphs.pkl)",
)
@click.option(
    "-o",
    "--outdir",
    default="pathways_plots",
    help="Output directory for plots",
    show_default=True,
)
@click.option(
    "--use-pydot",
    is_flag=True,
    help="Use pydot instead of graphviz",
)
@click.version_option(version=get_version(), prog_name="plot_modules_graphs")
def main(
    input_completeness,
    input_modules_list,
    modules_file,
    file_separator,
    graphs,
    outdir,
    use_pydot,
):
    """
    Generate plots for KEGG module pathways.

    Can visualize modules with or without completeness information.
    Requires at least one of: --input-completeness, --modules, or --modules-file.

    Examples:
    \b
    # Plot with completeness data
    plot_modules_graphs -i completeness.tsv -g graphs.pkl

    \b
    # Plot specific modules
    plot_modules_graphs -m M00001 -m M00002 -g graphs.pkl

    \b
    # Plot modules from file
    plot_modules_graphs -l modules.txt -g graphs.pkl
    """
    if not input_completeness and not input_modules_list and not modules_file:
        raise click.UsageError(
            "Must provide at least one of: --input-completeness, --modules, or --modules-file"
        )

    # Handle mutually exclusive options for modules input
    if input_modules_list and modules_file:
        raise click.UsageError("Cannot use both --modules and --modules-file")

    # Get graphs file
    graphs_filename = (
        graphs
        if graphs
        else files("kegg_pathways_completeness.pathways_data").joinpath("graphs.pkl")
    )

    plot_completeness_generator = PlotModuleCompletenessGraph(
        modules_completeness=parse_completeness_input(input_completeness),
        graphs=parse_graphs_input(graphs_filename),
        modules_definitions=None,
        outdir=outdir,
        modules_list=parse_input_modules(
            list(input_modules_list) if input_modules_list else None,
            modules_file,
            file_separator,
        ),
        use_pydot=use_pydot,
    )
    plot_completeness_generator.generate_plot()


if __name__ == "__main__":
    main()
