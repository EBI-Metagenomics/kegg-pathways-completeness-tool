#!/usr/bin/env python3

import argparse
import sys
import pickle
import networkx as nx
import logging
import copy
import os
from importlib.resources import files

from .plot_completeness_graphs import PlotModuleCompletenessGraph

__version__ = "1.1.0"

def setup_logging(verbose):
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s %(levelname)s - %(message)s'
    )


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Script generates modules pathways completeness by given set of KOs")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    input_args = parser.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-i", "--input", dest="input_file", help="Each line = pathway")
    input_args.add_argument("-l", "--input-list", dest="input_list", help="File with KOs comma separated")
    input_args.add_argument("-s", "--list-separator", dest="list_separator",
                            help="Separator for list option (default: comma)", default=',')

    parser.add_argument("-g", "--graphs", dest="graphs", help="graphs in pickle format", required=False)
    parser.add_argument("-a", "--pathways", dest="pathways", help="Pathways list", required=False)
    parser.add_argument("-n", "--names", dest="names", help="Pathway names", required=False)
    parser.add_argument("-c", "--classes", dest="classes", help="Pathway classes", required=False)

    parser.add_argument("-o", "--outdir", dest="outdir", help="output directory", default=".")
    parser.add_argument("-r", "--outprefix", dest="outprefix", help="prefix for output filename", default="summary.kegg")
    parser.add_argument("-w", "--include-weights", dest="include_weights", help="add weights for each KO in output",
                        action='store_true')
    parser.add_argument("-p", "--plot-pathways", dest="plot_pathways", help="Create images with pathways completeness",
                        action='store_true')
    parser.add_argument("-v", "--verbose", dest="verbose", help="Print more logging", required=False,
                        action='store_true')
    return parser.parse_args(argv)


def intersection(lst1, lst2):
    """
    Intersection of two sets
    :param lst1: first input list
    :param lst2: second input list
    :return: intersection in list format
    """
    return list(set(lst1) & set(lst2))


class CompletenessCalculator():
    def __init__(
            self,
            input_table: str,
            input_list: str,
            list_separator: str,
            outdir: str,
            outprefix: str,
            include_weights: bool,
            plot_pathways: bool,
            graphs: str = None,
            modules_definitions: str = None,
            module_classes: str = None,
            module_names: str = None,
    ):
        self.input_table = input_table
        self.input_list = input_list
        self.list_separator = list_separator

        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.path.makedirs(self.outdir)
        self.outprefix = outprefix
        self.name_output = os.path.join(self.outdir, self.outprefix)
        self.name_common_output_summary = os.path.join(self.name_output + '_pathways.tsv')
        self.name_contigs_output_summary = os.path.join(self.name_output + '_contigs.tsv')

        self.include_weights = include_weights
        self.plot_pathways = plot_pathways

        # load graphs
        self.graphs_name = graphs if graphs else files('kegg_pathways_completeness.graphs').joinpath('graphs.pkl')
        if os.path.exists(self.graphs_name):
            with open(self.graphs_name, 'rb') as graph_file:
                self.graphs = pickle.load(graph_file)
        else:
            logging.error(f'No graphs file found in {self.graphs_name}')
            sys.exit(1)

        # load pathways definitions
        self.modules_definitions_name = modules_definitions if modules_definitions \
            else files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways.txt')
        self.modules_definitions = {}
        if os.path.exists(self.modules_definitions_name):
            with open(self.modules_definitions_name, 'r') as pathways_file:
                for line in pathways_file:
                    fields = line.strip().split(':')
                    self.modules_definitions[fields[0]] = fields[1]
        else:
            logging.error(f'No {self.modules_definitions_name} found')
            sys.exit(1)

        # parse classes
        module_classes_name = module_classes if module_classes \
            else files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways_class.txt')
        self.module_classes = {}
        if os.path.exists(module_classes_name):
            with open(module_classes_name, 'r') as file_names:
                for line in file_names:
                    line = line.strip().split(':')
                    self.module_classes[line[0]] = line[1]
        else:
            logging.error(f'No {module_classes_name} found')
            sys.exit(1)

        # parse names
        module_names_name = module_names if module_names \
            else files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways_names.txt')
        self.module_names = {}
        if os.path.exists(module_names_name):
            with open(module_classes_name, 'r') as file_classes:
                for line in file_classes:
                    line = line.strip().split(':')
                    self.module_names[line[0]] = line[1]
        else:
            logging.error(f'No {module_names_name} found')
            sys.exit(1)

        self.edges, self.dict_KO_by_contigs = self.get_list_items()
        # !!!! careful - it was used a deepcopy
        self.weights_of_KOs = self.get_weights_for_KOs(self.graphs)

    def get_list_items(self):
        """
        Function creates a list of items that were found by HMMScan
        :param input_path: file with contigs and their KEGG annotations
        :return: list of unique KOs, { contig_name1: [KO1, KO2,...], contig_name2: [...], ...}
        """
        if not self.input_table and not self.input_list:
            logging.error("No necessary input table provided")
            sys.exit(1)
        items = []
        dict_KO_by_contigs = {}
        if self.input_table:
            if os.path.exists(self.input_table):
                with open(self.input_table, 'r') as file_in:
                    for line in file_in:
                        line = line.strip().split('\t')
                        name = line[0]
                        if name not in dict_KO_by_contigs:
                            dict_KO_by_contigs[name] = []
                        dict_KO_by_contigs[name] += line[1:]
                        items += line[1:]
            else:
                logging.error(f"No file {self.input_table}")
                sys.exit(1)
        elif self.input_list:
            if os.path.exists(self.input_list):
                name = os.path.basename(self.input_list)
                with open(self.input_list, 'r') as f:
                    list_kos = f.read().strip().split(self.list_separator)
                    if len(list_kos) == 0:
                        logging.error(f"No KOs found in {input_list}")
                    else:
                        dict_KO_by_contigs[name] = list_kos
                        items = list_kos
            else:
                logging.error(f"No file {self.input_list}")
                sys.exit(1)
        else:
            logging.error("No KOs provided")
        logging.info('KOs loaded')
        return list(set(items)), dict_KO_by_contigs

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
        dict_nodes_paths, dict_of_paths_labels, dict_of_weights, dict_of_new_weights = [{} for _ in range(4)]
        sorted_nodes = list(nx.topological_sort(G))
        for node in sorted_nodes:
            number_of_records = 0
            dict_nodes_paths[node], dict_of_paths_labels[node], dict_of_weights[node], dict_of_new_weights[node] \
                = [[], {}, {}, {}]
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
                    """ 
                        for multi edge pred---A---->node
                                         \____B____/ 
                    """
                    for num in range(number_of_pred_ancestors):
                        cur_labels = dict_of_paths_labels[pred][num]
                        dict_of_paths_labels[node][number_of_records] = \
                            cur_labels + [preds[pred][ancestor]['label']]
                        dict_of_weights[node][number_of_records] = \
                            dict_of_weights[pred][num] + preds[pred][ancestor]['weight']
                        dict_of_new_weights[node][number_of_records] = \
                            dict_of_new_weights[pred][num] + preds[pred][ancestor]['weight_new']
                        number_of_records += 1
                    for cur_path in dict_nodes_paths[pred]:
                        new_path = cur_path+[node]
                        dict_nodes_paths[node].append(new_path)
        paths_nodes, paths_labels = [dict_nodes_paths[1], dict_of_paths_labels[1]]
        weights, new_weights = [dict_of_weights[1], dict_of_new_weights[1]]
        metrics = []
        for num in range(len(weights)):
            metrics.append(1. * new_weights[num]/weights[num])
        indexes_min = [index for index in range(len(metrics)) if metrics[index] == min(metrics)]
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
                            if graph[start][finish][num]['label'] == edge:
                                index = num
                                break
                    else:
                        index = 0
                    # if graph[start][finish][index]['weight'] == 0:  # UNnecessary node
                    #     graph[start][finish][index]['weight'] = 1
                    graph[start][finish][index]['weight_new'] = 0

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
        percentage = round((1 - 1. * metrics[num]) * 100, 2)
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
        return percentage, len(indexes_min), list(matching_set), list(missing_set_necessary)


    def sort_out_pathways(self, contig_name, file_out_summary, edges):
        """
        Function sorts out all pathways and prints info about pathway that percentage of intersection more than 0
        :param
        contig_name == name of contig, or '' for full summary
        file_out_summary: output file
        :return: -
        """
        dict_sort_by_percentage = {}
        for name_pathway in self.graphs:
            graph = self.graphs[name_pathway]
            if intersection(graph[1], edges) == []:
                continue
            else:
                percentage, number_paths, matching_labels, missing_labels = \
                    self.calculate_percentage(graph=graph[0], dict_edges=graph[1], unnecessary_nodes=graph[2],
                                              edges=edges)
                if percentage != None:
                    if percentage not in dict_sort_by_percentage:
                        dict_sort_by_percentage[percentage] = {}
                    dict_sort_by_percentage[percentage][name_pathway] = [number_paths, matching_labels, missing_labels]

        # output Summary
        for percentage in sorted(list(dict_sort_by_percentage.keys()), reverse=True):
            #file_out_summary.write('**********************************************\nPercentage = ' + str(percentage) + '\n')
            for name_pathway in dict_sort_by_percentage[percentage]:
                matching_list = sorted(dict_sort_by_percentage[percentage][name_pathway][1])
                missing_list = sorted(dict_sort_by_percentage[percentage][name_pathway][2])
                if self.include_weights:
                    # matching
                    out_str = []
                    for KO in matching_list:
                        record = KO + '(' + str(weights_of_KOs[name_pathway][KO]) + ')'
                        out_str.append(record)
                    matching_current = ','.join(out_str)
                    # missing
                    out_str = []
                    for KO in missing_list:
                        out_str.append(KO + '(' + str(weights_of_KOs[name_pathway][KO]) + ')')
                    missing_current = ','.join(out_str)
                else:
                    matching_current = ','.join(matching_list)
                    missing_current = ','.join(missing_list)

                if contig_name != '':
                    out_name_pathway = '\t'.join([contig_name, name_pathway])
                else:
                    out_name_pathway = name_pathway
                output_line = '\t'.join([out_name_pathway, str(percentage), self.module_names[name_pathway],
                                        self.module_classes[name_pathway], matching_current, missing_current])
                file_out_summary.write(output_line + '\n')
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
        summary_header = '\t'.join(['module_accession', 'completeness', 'pathway_name',
                                              'pathway_class', 'matching_ko', 'missing_ko'])
        if contig:
            summary_header = 'contig\t' + summary_header
        file_summary.write(summary_header + '\n')

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
                        KO = graph[0]._adj[start][finish][num]['label']
                        weight = round(graph[0]._adj[start][finish][num]['weight'], 2)
                        dict_graphKO[name_pathway][KO] = weight
        logging.info('weights done')
        return dict_graphKO

    def generate_common_summary(self):
        # COMMON INFO
        logger = logging.getLogger(__name__)
        logger.info('Generating completeness for whole list of KOs...')
        #using_graphs = copy.deepcopy(graphs)
        with open(self.name_common_output_summary, "wt") as file_out_summary:
            self.set_headers(file_out_summary, contig=False)
            self.sort_out_pathways(contig_name='', file_out_summary=file_out_summary, edges=self.edges)
        logger.info('...Done')

    def generate_per_contig_summary(self):
        logger = logging.getLogger(__name__)
        logger.info('Generating completeness for contigs...')
        with open(self.name_contigs_output_summary, "wt") as file_out_summary:
            self.set_headers(file_out_summary, contig=True)
            for contig in self.dict_KO_by_contigs:
                #using_graphs = copy.deepcopy(graphs)
                edges = self.dict_KO_by_contigs[contig]
                self.sort_out_pathways(contig_name=contig, file_out_summary=file_out_summary, edges=edges)
        logger.info('...Done')

    def process(self):
        logger = logging.getLogger(__name__)
        # summary for all contigs
        self.generate_common_summary()
        # plot
        if self.plot_pathways:
            logger.info('Plot pathways images')
            plot_completeness_generator = PlotModuleCompletenessGraph(
                completeness_file=self.name_common_output_summary,
                graphs_pkl=self.graphs_name,
                modules_list=self.modules_definitions_name
            )
            plot_completeness_generator.generate_plot_for_completeness()
            logger.info('...Done. Results are in pathways_plots folder')

        # generate summary per-contig
        self.generate_per_contig_summary()
        logger.info('Bye!')


def main():
    args = parse_args(sys.argv[1:])
    setup_logging(args.verbose)

    completeness_calculator = CompletenessCalculator(
        input_table=args.input_file,
        input_list=args.input_list,
        list_separator=args.list_separator,
        outdir=args.outdir,
        outprefix=args.outprefix,
        include_weights=args.include_weights,
        plot_pathways=args.plot_pathways,
        graphs=args.graphs,
        module_names=args.names,
        module_classes=args.classes,
        modules_definitions=args.pathways
    )

    completeness_calculator.process()


if __name__ == "__main__":
    main()
