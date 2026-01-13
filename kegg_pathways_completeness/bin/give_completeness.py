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

import sys
import networkx as nx
import logging
import copy
import os
import click
from importlib.resources import files
from importlib.metadata import version, PackageNotFoundError

from .plot_modules_graphs import PlotModuleCompletenessGraph
from .utils import (
    parse_graphs_input,
    setup_logging,
    intersection,
)


def get_version():
    """Get package version from installed metadata"""
    try:
        return version("kegg-pathways-completeness")
    except PackageNotFoundError:
        return "unknown"


def parse_modules_table_tsv(tsv_file):
    """
    Parse new TSV format file (modules_table.tsv) into separate dictionaries.
    Returns: (modules_definitions, modules_names, modules_classes)
    """
    modules_definitions = {}
    modules_names = {}
    modules_classes = {}

    with open(tsv_file, "r") as f:
        # Read header
        header = f.readline().strip().split("\t")

        # Validate header
        required_cols = ["module", "definition", "name", "class"]
        for col in required_cols:
            if col not in header:
                raise ValueError(f"TSV file must have '{col}' column")

        module_idx = header.index("module")
        definition_idx = header.index("definition")
        name_idx = header.index("name")
        class_idx = header.index("class")

        # Read data rows
        for line in f:
            if line.strip():
                fields = line.strip().split("\t")
                if len(fields) > max(module_idx, definition_idx, name_idx, class_idx):
                    module = fields[module_idx]
                    modules_definitions[module] = fields[definition_idx]
                    modules_names[module] = fields[name_idx]
                    modules_classes[module] = fields[class_idx]

    return modules_definitions, modules_names, modules_classes


class CompletenessCalculator:
    def __init__(
        self,
        input_KOs: dict,
        outdir: str,
        outprefix: str,
        include_weights: bool,
        plot_pathways: bool,
        per_contig: bool,
        graphs: nx.MultiDiGraph,
        modules_definitions: str = None,
        modules_classes: str = None,
        modules_names: str = None,
    ):
        # input KOs
        self.dict_KO_by_contigs = input_KOs

        # modules graphs
        self.graphs = graphs

        # modules info
        self.modules_definitions = modules_definitions
        self.modules_classes = modules_classes
        self.modules_names = modules_names

        # outputs
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.outprefix = outprefix
        self.name_output = os.path.join(self.outdir, self.outprefix)
        self.name_output_pathways_plots = os.path.join(self.outdir, "pathways_plots")

        # flags
        self.include_weights = include_weights
        if self.include_weights:
            self.name_common_output_summary = os.path.join(
                self.name_output + "_pathways.with_weights.tsv"
            )
            self.name_contigs_output_summary = os.path.join(
                self.name_output + "_contigs.with_weights.tsv"
            )
        else:
            self.name_common_output_summary = os.path.join(
                self.name_output + "_pathways.tsv"
            )
            self.name_contigs_output_summary = os.path.join(
                self.name_output + "_contigs.tsv"
            )
        self.plot_pathways = plot_pathways
        self.per_contig = per_contig

        self.edges = self.get_edges_list()
        # !!!! careful - it was used a deepcopy
        self.weights_of_KOs = self.get_weights_for_KOs(copy.deepcopy(self.graphs))

    def get_edges_list(self):
        items = []
        for contig in self.dict_KO_by_contigs:
            items += self.dict_KO_by_contigs[contig]
        return list(set(items))

    def finding_paths(self, G):
        """
        Function sorts out all paths in the given graph. Moreover, for each found path calculating metrics M.
        M = weight_new_path / weigth_old_path --> min
        :param G: graph
        :return: paths_nodes - sequence of nodes that represents path (ex. [0 2 4 5 6 7 1])
                 paths_labels - sequence of labels of items that represents path (ex. [K1 K3 K4 K0])
                 weights - old weights of paths
                 new_weights - new weights of paths
                 indexes_min - list of indexes of paths with the smallest M
        """
        dict_nodes_paths, dict_of_paths_labels, dict_of_weights, dict_of_new_weights = [
            {} for _ in range(4)
        ]
        sorted_nodes = list(nx.topological_sort(G))
        for node in sorted_nodes:
            number_of_records = 0
            (
                dict_nodes_paths[node],
                dict_of_paths_labels[node],
                dict_of_weights[node],
                dict_of_new_weights[node],
            ) = [[], {}, {}, {}]
            preds = G.pred[node]
            if preds == {}:
                dict_nodes_paths[node].append([node])
                dict_of_paths_labels[node][0] = []
                dict_of_weights[node][0] = 0
                dict_of_new_weights[node][0] = 0
                continue
            for pred in preds.keys():  # ancestors of node: pred-->node
                number_of_pred_ancestors = len(dict_of_paths_labels[pred])

                for ancestor in preds[pred]:
                    # for multi edge pred---A---->node
                    #                      \____B____/
                    for num in range(number_of_pred_ancestors):
                        cur_labels = dict_of_paths_labels[pred][num]
                        dict_of_paths_labels[node][number_of_records] = cur_labels + [
                            preds[pred][ancestor]["label"]
                        ]
                        dict_of_weights[node][number_of_records] = (
                            dict_of_weights[pred][num] + preds[pred][ancestor]["weight"]
                        )
                        dict_of_new_weights[node][number_of_records] = (
                            dict_of_new_weights[pred][num]
                            + preds[pred][ancestor]["weight_new"]
                        )
                        number_of_records += 1
                    for cur_path in dict_nodes_paths[pred]:
                        new_path = cur_path + [node]
                        dict_nodes_paths[node].append(new_path)
        paths_nodes, paths_labels = [dict_nodes_paths[1], dict_of_paths_labels[1]]
        weights, new_weights = [dict_of_weights[1], dict_of_new_weights[1]]
        metrics = []
        for num in range(len(weights)):
            metrics.append(1.0 * new_weights[num] / weights[num])
        indexes_min = [
            index for index in range(len(metrics)) if metrics[index] == min(metrics)
        ]
        return paths_nodes, paths_labels, metrics, indexes_min

    def calculate_percentage(self, graph, dict_edges, unnecessary_nodes, edges):
        """
        Function returns the percentage of matches of set of edges and graph.
        Example:
                Pathway: A B C. Edges: A -> percentage = 33
        :param graph: input graph of pathway
        :param dict_edges: dict of edges in graph by labels
        :param unnecessary_nodes: list of all nodes
        :return: percentage [0:100], number of paths, matching_list, missing_list of KOs
        """
        # set weights_new as 0 for edges that are presented
        for edge in edges:
            if edge in dict_edges:
                nodes = dict_edges[edge]
                for cur_pair in nodes:
                    start = cur_pair[0]
                    finish = cur_pair[1]
                    if len(graph[start][finish]) > 0:
                        for num in range(len(graph[start][finish])):
                            if graph[start][finish][num]["label"] == edge:
                                index = num
                                break
                    else:
                        index = 0
                    # if graph[start][finish][index]['weight'] == 0:  # UNnecessary node
                    #     graph[start][finish][index]['weight'] = 1
                    graph[start][finish][index]["weight_new"] = 0

        # find the best path(s)
        paths_nodes, paths_labels, metrics, indexes_min = self.finding_paths(graph)
        """
        metrics [N]: list of all sum_weight_new/sum_weight for all possible paths in graph
        paths_nodes [N]: list of nodes that construct each path
        paths_labels [N]: list of labels that construct each path
        indexes_min [M<=N]: list of indexes that correspond to minimum value in metrict
        """

        # take random path from minimal, for example first
        # because all paths in indexes_min have the same percentage. That means there is no difference which one to output
        num = indexes_min[0]
        percentage = round((1 - 1.0 * metrics[num]) * 100, 2)
        matching_set, missing_set_necessary, missing_set = [set() for _ in range(3)]
        if percentage > 0:
            new_labels = paths_labels[num]
            missing_labels = set(new_labels).difference(set(edges))
            missing_set = missing_set.union(missing_labels)
            missing_set_necessary = missing_set.difference(set(unnecessary_nodes))

            existing_labels = set(new_labels).intersection(set(edges))
            matching_set = matching_set.union(existing_labels)
        else:
            percentage = None
        return (
            percentage,
            len(indexes_min),
            list(matching_set),
            list(missing_set_necessary),
        )

    def sort_out_pathways(self, contig_name, file_out_summary, edges, graphs):
        """
        Function sorts out all pathways and prints info about pathway that percentage of intersection more than 0
        :param
        contig_name == name of contig, or '' for full summary
        file_out_summary: output file
        :return: -
        """
        dict_sort_by_percentage, module_matching_kos = {}, {}
        for name_pathway in graphs:
            graph = graphs[name_pathway]
            if intersection(graph[1], edges) == []:
                continue
            else:
                (
                    percentage,
                    number_paths,
                    matching_labels,
                    missing_labels,
                ) = self.calculate_percentage(
                    graph=graph[0],
                    dict_edges=graph[1],
                    unnecessary_nodes=graph[2],
                    edges=edges,
                )
                if percentage is not None:
                    if percentage not in dict_sort_by_percentage:
                        dict_sort_by_percentage[percentage] = {}
                    dict_sort_by_percentage[percentage][name_pathway] = [
                        number_paths,
                        matching_labels,
                        missing_labels,
                    ]

        # output Summary
        for percentage in sorted(list(dict_sort_by_percentage.keys()), reverse=True):
            # file_out_summary.write('**********************************************\nPercentage = ' + str(percentage) + '\n')
            for name_pathway in dict_sort_by_percentage[percentage]:
                matching_list = sorted(
                    dict_sort_by_percentage[percentage][name_pathway][1]
                )
                missing_list = sorted(
                    dict_sort_by_percentage[percentage][name_pathway][2]
                )
                if self.include_weights:
                    # matching
                    out_str = []
                    for KO in matching_list:
                        record = (
                            KO + "(" + str(self.weights_of_KOs[name_pathway][KO]) + ")"
                        )
                        out_str.append(record)
                    matching_current = ",".join(out_str)
                    # missing
                    out_str = []
                    for KO in missing_list:
                        out_str.append(
                            KO + "(" + str(self.weights_of_KOs[name_pathway][KO]) + ")"
                        )
                    missing_current = ",".join(out_str)
                else:
                    matching_current = ",".join(matching_list)
                    missing_current = ",".join(missing_list)

                if contig_name != "":
                    out_name_pathway = "\t".join([contig_name, name_pathway])
                else:
                    out_name_pathway = name_pathway
                module_matching_kos[name_pathway] = matching_current
                output_line = "\t".join(
                    [
                        out_name_pathway,
                        str(percentage),
                        self.modules_names[name_pathway],
                        self.modules_classes[name_pathway],
                        matching_current,
                        missing_current,
                    ]
                )
                file_out_summary.write(output_line + "\n")
        return module_matching_kos
        """
        file_out_summary.write('\n******* REMINDER ********')
        file_out_summary.write('Number of nodes: ' + str(len(edges)) + '\n')
        file_out_summary.write('Set of nodes: ' + str(edges) + '\n')
        """

    def set_headers(self, file_summary, contig):
        """
        Function adds header to output result tables.
        :param file_summary: output file
        :param contig: flag true (running for every contig) or false (running in general mode)
        :return: -
        """
        summary_header = "\t".join(
            [
                "module_accession",
                "completeness",
                "pathway_name",
                "pathway_class",
                "matching_ko",
                "missing_ko",
            ]
        )
        if contig:
            summary_header = "contig\t" + summary_header
        file_summary.write(summary_header + "\n")

    def get_weights_for_KOs(self, graphs):
        """
        For each graph functions returns dict of weights by KOs,
        ex. dict_graphKO: { pathway1: {KO1: weight1, KO2: weight2}, pathway2: {...} }
        :param graphs: dict of graphs
        :return: dict of pathways with weights for each KO
        """
        dict_graphKO = {}
        for name_pathway in graphs:
            graph = graphs[name_pathway]
            dict_graphKO[name_pathway] = {}

            for start in graph[0]._adj:
                for finish in graph[0]._adj[start]:
                    for num in graph[0]._adj[start][finish]:
                        KO = graph[0]._adj[start][finish][num]["label"]
                        weight = round(graph[0]._adj[start][finish][num]["weight"], 2)
                        dict_graphKO[name_pathway][KO] = weight
        logging.info("weights done")
        return dict_graphKO

    def generate_common_summary(self):
        # COMMON INFO
        logger = logging.getLogger(__name__)
        logger.info("Generating completeness for whole list of KOs...")
        using_graphs = copy.deepcopy(self.graphs)
        with open(self.name_common_output_summary, "wt") as file_out_summary:
            self.set_headers(file_out_summary, contig=False)
            module_matching_kos = self.sort_out_pathways(
                contig_name="",
                file_out_summary=file_out_summary,
                edges=self.edges,
                graphs=using_graphs,
            )
        logger.info("...Done")
        return module_matching_kos

    def generate_per_contig_summary(self):
        logger = logging.getLogger(__name__)
        logger.info("Generating completeness for contigs...")
        with open(self.name_contigs_output_summary, "wt") as file_out_summary:
            self.set_headers(file_out_summary, contig=True)
            for contig in self.dict_KO_by_contigs:
                using_graphs = copy.deepcopy(self.graphs)
                edges = self.dict_KO_by_contigs[contig]
                self.sort_out_pathways(
                    contig_name=contig,
                    file_out_summary=file_out_summary,
                    edges=edges,
                    graphs=using_graphs,
                )
        logger.info("...Done")

    def process(self):
        logger = logging.getLogger(__name__)
        # summary for all contigs
        module_matching_kos = self.generate_common_summary()
        # plot
        if self.plot_pathways:
            logger.info("Plot pathways images")
            plot_completeness_generator = PlotModuleCompletenessGraph(
                modules_completeness=module_matching_kos,
                graphs=self.graphs,
                modules_definitions=self.modules_definitions,
                outdir=self.name_output_pathways_plots,
            )
            plot_completeness_generator.generate_plot()
            logger.info("...Done. Results are in pathways_plots folder")

        # generate summary per-contig
        if self.per_contig:
            self.generate_per_contig_summary()
        logger.info("Bye!")


def get_kos_dict(input_table, input_list, list_separator):
    """
    Function creates a list of items that were found by HMMScan
    :param input_path: file with contigs and their KEGG annotations
    :return: list of unique KOs, { contig_name1: [KO1, KO2,...], contig_name2: [...], ...}
    """
    if not input_table and not input_list:
        logging.error("No necessary input table provided")
        sys.exit(1)
    dict_KO_by_contigs = {}
    if input_table:
        if os.path.exists(input_table):
            with open(input_table, "r") as file_in:
                for line in file_in:
                    line = line.strip().split("\t")
                    name = line[0]
                    if name not in dict_KO_by_contigs:
                        dict_KO_by_contigs[name] = []
                    dict_KO_by_contigs[name] += line[1:]
        else:
            logging.error(f"No file {input_table}")
            sys.exit(1)
    elif input_list:
        if os.path.exists(input_list):
            name = os.path.basename(input_list)
            with open(input_list, "r") as f:
                list_kos = f.read().strip().split(list_separator)
                if len(list_kos) == 0:
                    logging.error(f"No KOs found in {input_list}")
                else:
                    dict_KO_by_contigs[name] = list_kos
        else:
            logging.error(f"No file {input_list}")
            sys.exit(1)
    else:
        logging.error("No KOs provided")
    logging.info("KOs loaded")
    return dict_KO_by_contigs


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True),
    help="Input file with KOs (each line = contig_name\\tKO1\\tKO2\\t...)",
)
@click.option(
    "-l",
    "--input-list",
    "input_list",
    type=click.Path(exists=True),
    help="File with KOs separated by delimiter (see --list-separator)",
)
@click.option(
    "-s",
    "--list-separator",
    default=",",
    help="Separator for --input-list option",
    show_default=True,
)
@click.option(
    "-g",
    "--graphs",
    type=click.Path(exists=True),
    help="Graphs in pickle format (default: uses packaged graphs.pkl)",
)
@click.option(
    "-t",
    "--modules-table",
    type=click.Path(exists=True),
    help="Modules table in TSV format (default: uses packaged modules_table.tsv)",
)
@click.option(
    "-o",
    "--outdir",
    default=".",
    help="Output directory",
    show_default=True,
)
@click.option(
    "-r",
    "--outprefix",
    default="summary.kegg",
    help="Prefix for output filename",
    show_default=True,
)
@click.option(
    "-w",
    "--include-weights",
    is_flag=True,
    help="Add weights for each KO in output",
)
@click.option(
    "-p",
    "--plot-pathways",
    is_flag=True,
    help="Create images with pathways completeness",
)
@click.option(
    "-m",
    "--add-per-contig",
    "per_contig",
    is_flag=True,
    help="Create per-contig summary table (does not work with --input-list)",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Enable verbose logging",
)
@click.version_option(version=get_version(), prog_name="give_completeness")
def main(
    input_file,
    input_list,
    list_separator,
    graphs,
    modules_table,
    outdir,
    outprefix,
    include_weights,
    plot_pathways,
    per_contig,
    verbose,
):
    """
    Calculate KEGG pathway completeness for a given set of KOs.

    Input options (mutually exclusive):
    \b
    - Use -i/--input for a file with contigs and KOs
    - Use -l/--input-list for a simple list of KOs

    Module data:
    \b
    - Use -t/--modules-table for modules data in TSV format (required)
    """
    setup_logging(verbose)
    logger = logging.getLogger(__name__)

    # Validate mutually exclusive inputs
    if not input_file and not input_list:
        raise click.UsageError("Must provide either --input or --input-list")
    if input_file and input_list:
        raise click.UsageError("Cannot use both --input and --input-list")

    # Parse input with KOs
    dict_KO_by_contigs = get_kos_dict(
        input_table=input_file,
        input_list=input_list,
        list_separator=list_separator,
    )

    # Get graphs file
    graphs_filename = (
        graphs
        if graphs
        else files("kegg_pathways_completeness.pathways_data").joinpath("graphs.pkl")
    )

    # Get modules table file
    modules_table_filename = (
        modules_table
        if modules_table
        else files("kegg_pathways_completeness.pathways_data").joinpath("modules_table.tsv")
    )

    # Load modules information from TSV
    logger.info(f"Loading modules data from TSV: {modules_table_filename}")
    modules_definitions, modules_names, modules_classes = parse_modules_table_tsv(
        modules_table_filename
    )

    completeness_calculator = CompletenessCalculator(
        input_KOs=dict_KO_by_contigs,
        outdir=outdir,
        outprefix=outprefix,
        graphs=parse_graphs_input(graphs_filename),
        modules_names=modules_names,
        modules_classes=modules_classes,
        modules_definitions=modules_definitions,
        include_weights=include_weights,
        plot_pathways=plot_pathways,
        per_contig=per_contig,
    )

    completeness_calculator.process()


if __name__ == "__main__":
    main()
