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
import time

import click
import networkx as nx
import numpy as np

from .utils import get_version, setup_logging


class GraphsGenerator:
    def __init__(
        self,
        input_file: str,
        output_dir: str,
        existing_graphs_file: str = None,
        changed_file: str = None,
    ):
        """
        Creates Graphs in network format for each module in input_file.
        Graphs object is saved into graphs.pkl file.

        Can perform incremental updates by reusing existing graphs and only
        regenerating changed modules.

        :param input_file: Line separated file with modules in format module:KOs
        :param output_dir: name of output directory
        :param existing_graphs_file: Path to existing graphs.pkl to reuse unchanged graphs
        :param changed_file: Path to changed.tsv with modules that need regeneration
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.existing_graphs_file = existing_graphs_file
        self.changed_file = changed_file
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

    def set_order_separators(self, dict_levels):
        """
        Order of parsing: [',', ' ', '+', '-']
        Function returns the parsing separators in the right order according to levels
        ??? not optimal ???
        """
        keys = sorted(list(dict_levels.keys()))
        if keys != []:
            min_level, max_level = [
                int(keys[0].split("_")[0]),
                int(keys[len(keys) - 1].split("_")[0]) + 1,
            ]
            orders = [
                str(j) + "_" + i
                for j in range(min_level, max_level)
                for i in [",", " ", "+", "-"]
            ]
            new_order = [element for element in orders if element in keys]
        else:
            new_order = []
        return new_order

    def add_to_dict_of_levels(self, dict_levels, c, cur_level, index):
        """
        Function returns the dict of positions according to the level of space or comma
        Example: {'1_,': [14], '2_,': [9], '0_ ': [3], '1_ ': [12]}
            comma of level 1: position 14
            comma of level 2: position 9
            space of level 0: position 3
            space of level 1: position 12
        """
        symbol = str(cur_level) + "_" + c
        if symbol not in dict_levels:
            dict_levels[symbol] = []
        dict_levels[symbol].append(index)
        return dict_levels

    def set_brackets(self, pathway):
        """
        Function defines levels of all brackets in expression. The output will be used by function <check_brackets>
        Example 1:
            expression: A B (C,D)
            levels:  -1,-1,-1,-1,0,-1,-1,-1,0
        Example 2:
            expression: (A B (C,D))
            levels:  0,-1,-1,-1,-1,1,-1,-1,-1,1,0
        :param pathway: string expression
        :return: levels of brackets
        """
        levels_brackets = []
        cur_open = []
        num = -1
        for c in pathway:
            if c == "(":
                num += 1
                cur_open.append(num)
                levels_brackets.append(num)
            elif c == ")":
                levels_brackets.append(cur_open[len(cur_open) - 1])
                cur_open.pop()
            else:
                levels_brackets.append(-1)
        return levels_brackets

    def set_levels(self, pathway):
        """
        Function creates a dictionary of separators in pathway.
           Keys format: level_separator (ex. '1_,' or '0_ ')
           Values: list of positions in expression
        Example:
            expression: D (A+B) -> levels: 0011111 -> dict_levels: {'0_ ':[1], '1+':[4] }

        :param pathway: string expression
        :return: dict. of separators with their positions
        """
        dict_levels = {}
        L = len(pathway)
        cur_level, index = [0 for _ in range(2)]

        while index < L:
            c = pathway[index]
            if c == " " or c == "," or c == "-" or c == "+":
                dict_levels = self.add_to_dict_of_levels(
                    dict_levels, c, cur_level, index
                )
            elif c == "(":
                cur_level += 1
            elif c == ")":
                cur_level -= 1
            else:
                index += 1
                if index < L:
                    while pathway[index] not in [" ", ",", "(", ")", "-", "+"]:
                        index += 1
                        if index >= L:
                            break
                    index -= 1
            index += 1
        return dict_levels

    def check_brackets(self, pathway, levels_brackets):
        """
        Function checks is this expression in brackets. Returns without if true
        Example: input (A B C) or ((A B C))
                return: A B C
        :param pathway: input string expression
        :return: output string expression
        """
        expression_without_brackets = pathway
        L = len(pathway)
        for i in range(L):
            # check brackets
            if (
                pathway[i] == "("
                and pathway[L - i - 1] == ")"
                and levels_brackets[i] == levels_brackets[L - i - 1]
            ):
                expression_without_brackets = pathway[i + 1: L - i - 1]
            else:
                return expression_without_brackets
        return expression_without_brackets

    def recursive_parsing(
        self, G, dict_edges, unnecessary_nodes, expression, start_node, end_node, weight
    ):
        """
        Main parser:
           - adds edges and nodes to global graph G
           - adds names of edges to global dictionary of edges

        :param expression: current string expression to parse
        :param start_node: num of node from which expression sequence would be started
        :param end_node: num of node to which expression sequence would be finished
        :param weight: weight of edge (0 for unnecessary edges, 1 - for necessary, float - for parts of complex)
        :return: graph, dict of edges
        """
        if expression == "--":  # case --
            name_missing = "K00000"
            # print('MAKE EDGE: ' + name_missing, start_node, end_node)
            G.add_edge(
                start_node,
                end_node,
                label=name_missing,
                weight=0,
                weight_new=0,
                name="-",
            )
            unnecessary_nodes.append(name_missing)
            if name_missing not in dict_edges:
                dict_edges[name_missing] = []
            dict_edges[name_missing].append([start_node, end_node])
            return G, dict_edges, unnecessary_nodes

        expression = self.check_brackets(
            expression, self.set_brackets(expression)
        )  # delete brackets (expression)
        cur_dict_levels = self.set_levels(
            expression
        )  # define levels of each part of expression
        separators_order = self.set_order_separators(
            cur_dict_levels
        )  # set orders of existing separators
        cur_weight = weight

        if len(separators_order) == 1:  # case: -K....
            if separators_order[0] == "0_-" and expression[0] == "-":
                # print('MAKE EDGE: ' + expression[1:], start_node, end_node)
                G.add_edge(
                    start_node,
                    end_node,
                    label=expression[1:],
                    weight=0,
                    weight_new=0,
                    name="-",
                )
                unnecessary_nodes.append(expression[1:])
                if expression[1:] not in dict_edges:
                    dict_edges[expression[1:]] = []
                dict_edges[expression[1:]].append([start_node, end_node])
                return G, dict_edges, unnecessary_nodes

        if separators_order != []:
            # separator
            field = separators_order[0]
            symbol = field.split("_")[1]

            if symbol == "+" or symbol == " ":
                cur_weight = cur_weight / (len(cur_dict_levels[field]) + 1)

            separators = list(np.array(sorted(cur_dict_levels[field])))
            cur_sep = 0
            cur_start_node = start_node
            cur_end_node = end_node

            for separator, num in zip(separators, range(len(separators))):

                if symbol == " " or symbol == "+" or symbol == "-":
                    cur_end_node = len(list(G.nodes()))
                    G.add_node(cur_end_node)
                if symbol == "-" and num > 0:
                    cur_weight = 0

                subexpression = expression[cur_sep:separator]

                G, dict_edges, unnecessary_nodes = self.recursive_parsing(
                    G=G,
                    dict_edges=dict_edges,
                    unnecessary_nodes=unnecessary_nodes,
                    expression=subexpression,
                    start_node=cur_start_node,
                    end_node=cur_end_node,
                    weight=cur_weight,
                )
                cur_sep = separator + 1
                if symbol == " " or symbol == "+" or symbol == "-":
                    cur_start_node = cur_end_node

            num += 1
            if symbol == " " or symbol == "+" or symbol == "-":  # nodes and edges
                cur_start_node = cur_end_node
                cur_end_node = end_node
            if symbol == "-" and num > 0:  # weight
                cur_weight = 0

            G, dict_edges, unnecessary_nodes = self.recursive_parsing(
                G=G,
                dict_edges=dict_edges,
                unnecessary_nodes=unnecessary_nodes,
                expression=expression[cur_sep: len(expression)],
                start_node=cur_start_node,
                end_node=cur_end_node,
                weight=cur_weight,
            )
            return G, dict_edges, unnecessary_nodes
        else:
            # print('MAKE EDGE: ' + expression, start_node, end_node)
            if cur_weight == 0:
                G.add_edge(
                    start_node,
                    end_node,
                    label=expression,
                    weight=cur_weight,
                    weight_new=cur_weight,
                    name="-",
                )
                unnecessary_nodes.append(expression)
            else:
                G.add_edge(
                    start_node,
                    end_node,
                    label=expression,
                    weight=cur_weight,
                    weight_new=cur_weight,
                    name="node",
                )
            if expression not in dict_edges:
                dict_edges[expression] = []
            dict_edges[expression].append([start_node, end_node])
            return G, dict_edges, unnecessary_nodes

    def _is_tsv_format(self):
        """Detect if input file is in TSV format (new) or old format"""
        with open(self.input_file, "r") as f:
            first_line = f.readline().strip()
            # Check if first line is a TSV header
            return (
                first_line.startswith("module\t")
                or "\t" in first_line
                and ":" not in first_line
            )

    def _read_modules_from_tsv(self):
        """Read modules from new TSV format file"""
        modules = []
        with open(self.input_file, "r") as f:
            # Skip header
            header = f.readline().strip().split("\t")

            # Validate header
            if "module" not in header or "definition" not in header:
                raise ValueError("TSV file must have 'module' and 'definition' columns")

            module_idx = header.index("module")
            definition_idx = header.index("definition")

            # Read data rows
            for line in f:
                if line.strip():
                    fields = line.strip().split("\t")
                    if len(fields) > max(module_idx, definition_idx):
                        module = fields[module_idx]
                        definition = fields[definition_idx]
                        modules.append((module, definition))

        return modules

    def _read_modules_from_old_format(self):
        """Read modules from old format file (module:definition)"""
        modules = []
        with open(self.input_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and ":" in line:
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        modules.append((parts[0], parts[1]))
        return modules

    def _load_existing_graphs(self):
        """Load existing graphs from pickle file"""
        if not self.existing_graphs_file:
            return {}

        if not os.path.exists(self.existing_graphs_file):
            logging.warning(
                f"Existing graphs file not found: {self.existing_graphs_file}"
            )
            return {}

        try:
            with open(self.existing_graphs_file, "rb") as f:
                graphs = pickle.load(f)
            logging.info(
                f"Loaded {len(graphs)} existing graphs from {self.existing_graphs_file}"
            )
            return graphs
        except Exception as e:
            logging.error(f"Error loading existing graphs: {e}")
            return {}

    def _read_changed_modules(self):
        """Read list of changed modules from changed.tsv file"""
        if not self.changed_file:
            return set()

        if not os.path.exists(self.changed_file):
            logging.warning(f"Changed file not found: {self.changed_file}")
            return set()

        changed_modules = set()
        try:
            with open(self.changed_file, "r") as f:
                # Skip header
                f.readline()
                # Read module IDs from first column
                for line in f:
                    if line.strip():
                        fields = line.strip().split("\t")
                        if fields:
                            changed_modules.add(fields[0])
            logging.info(
                f"Found {len(changed_modules)} changed modules in {self.changed_file}"
            )
            return changed_modules
        except Exception as e:
            logging.error(f"Error reading changed modules: {e}")
            return set()

    def pathways_processing(self):
        """
        Main function for processing each pathway.
        Supports both old format (module:definition) and new TSV format.
        Supports incremental updates by reusing existing graphs.

        Function creates dictionary key: name; value: (graph, dict_edges)
        """
        logger = logging.getLogger(__name__)
        logger.info("Start graphs generation")

        # Detect file format and read modules
        if self._is_tsv_format():
            logger.info("Detected TSV format input file")
            modules = self._read_modules_from_tsv()
        else:
            logger.info("Detected old format input file (module:definition)")
            modules = self._read_modules_from_old_format()

        # Create dictionary for quick lookup
        modules_dict = {name: pathway for name, pathway in modules}
        logger.info(f"Total modules in input: {len(modules_dict)}")

        # Load existing graphs if provided
        existing_graphs = self._load_existing_graphs()

        # Read changed modules if provided
        changed_modules = self._read_changed_modules()

        # Determine which modules to process
        if existing_graphs and changed_modules:
            # Incremental mode: only regenerate changed modules
            modules_to_generate = changed_modules & set(modules_dict.keys())
            modules_to_reuse = set(modules_dict.keys()) - changed_modules

            logger.info("Incremental update mode:")
            logger.info(f"  - Modules to regenerate: {len(modules_to_generate)}")
            logger.info(f"  - Modules to reuse from existing: {len(modules_to_reuse)}")

            # Start with existing graphs for modules that haven't changed
            graphs = {}
            reused_count = 0
            for module in modules_to_reuse:
                if module in existing_graphs:
                    graphs[module] = existing_graphs[module]
                    reused_count += 1
                    logger.debug(f"Reusing graph for {module}")

            logger.info(f"Reused {reused_count} existing graphs")

            # Generate graphs only for changed modules
            modules_to_process = [
                (name, modules_dict[name]) for name in modules_to_generate
            ]
        else:
            # Full regeneration mode
            logger.info("Full regeneration mode (no incremental update)")
            graphs = {}
            modules_to_process = modules

        # Process modules (either all or just changed ones)
        if modules_to_process:
            logger.info(f"Generating {len(modules_to_process)} graphs...")
            for name, pathway in modules_to_process:
                logger.debug(f"Processing {name}")
                # Graph creation:
                Graph = nx.MultiDiGraph()
                Graph.add_node(0, color="green")
                Graph.add_node(1, color="red")
                # Parsing
                Graph, dict_edges, unnecessary_nodes = self.recursive_parsing(
                    G=Graph,
                    dict_edges={},
                    unnecessary_nodes=[],
                    expression=pathway,
                    start_node=0,
                    end_node=1,
                    weight=1,
                )
                # Saving
                graphs[name] = tuple([Graph, dict_edges, unnecessary_nodes])
                time.sleep(1)

        # Final summary
        logger.info(f"Total graphs in output: {len(graphs)}")
        if existing_graphs:
            removed = set(existing_graphs.keys()) - set(graphs.keys())
            if removed:
                logger.info(
                    f"Removed {len(removed)} graphs (no longer in input): {sorted(list(removed))[:10]}..."
                )

        # Save graphs
        logger.info("Done. Saving graphs...")
        path_output = os.path.join(self.output_dir, "graphs.pkl")
        with open(path_output, "wb") as f:
            pickle.dump(graphs, f)
        logger.info(f"Graphs saved to {path_output}")


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input file: TSV format (modules_table.tsv) or old format (module:definition per line)",
)
@click.option(
    "-o",
    "--outdir",
    default="outdir",
    help="Output directory where graphs.pkl will be stored",
    show_default=True,
)
@click.option(
    "-e",
    "--existing-graphs",
    type=click.Path(exists=True),
    help="Existing graphs.pkl file to reuse unchanged graphs (for incremental updates)",
)
@click.option(
    "-c",
    "--changed",
    type=click.Path(exists=True),
    help="Changed modules file (changed.tsv) - only regenerate graphs for these modules",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Enable verbose logging",
)
@click.version_option(version=get_version(), prog_name="make_graphs")
def main(input_file, outdir, existing_graphs, changed, verbose):
    """
    Generates graph structures for KEGG modules and saves them to graphs.pkl.

    Supports both formats:
    - New TSV format: modules_table.tsv with columns (module, definition, name, class)
    - Old format: one module per line as module:definition

    The script automatically detects the input format.

    Incremental Updates:
    \b
    Use --existing-graphs and --changed together for incremental updates:
    - Loads existing graphs from --existing-graphs
    - Only regenerates graphs for modules in --changed file
    - Reuses graphs for unchanged modules
    - Removes graphs for modules not in --input (deleted from KEGG)

    Examples:
    \b
    # Full regeneration
    make_graphs -i modules_table.tsv -o graphs_output

    \b
    # Incremental update (only changed modules)
    make_graphs -i modules_table.tsv -o graphs_output \\
                -e old_graphs.pkl -c changed.tsv
    """
    setup_logging(verbose)

    # Validate incremental update options
    if existing_graphs and not changed:
        logging.warning(
            "--existing-graphs provided without --changed. Will perform full regeneration."
        )
    if changed and not existing_graphs:
        logging.warning(
            "--changed provided without --existing-graphs. Will regenerate all modules in changed file."
        )

    graphs_generator = GraphsGenerator(
        input_file=input_file,
        output_dir=outdir,
        existing_graphs_file=existing_graphs,
        changed_file=changed,
    )
    graphs_generator.pathways_processing()


if __name__ == "__main__":
    main()
