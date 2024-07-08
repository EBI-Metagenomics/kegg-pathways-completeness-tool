#!/usr/bin/env python3

import argparse
import sys
import pickle
import networkx as nx
import logging
import copy
import os
from importlib.resources import files

from kegg_pathways_completeness.bin.plot_completeness_graphs import plot_graphs, parse_input

__version__ = "1.0.5"

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

def parse_args():
    default_graphs_path = files('kegg_pathways_completeness.graphs').joinpath('graphs.pkl')
    default_pathways_path = files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways.txt')
    default_pathways_names_path = files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways_names.txt')
    default_pathways_class_path = files('kegg_pathways_completeness.pathways_data').joinpath('all_pathways_class.txt')

    parser = argparse.ArgumentParser(description="Script generates modules pathways completeness by given set of KOs")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    input_args = parser.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-i", "--input", dest="input_file", help="Each line = pathway")
    input_args.add_argument("-l", "--input-list", dest="input_list", help="File with KOs comma separated")

    parser.add_argument("-g", "--graphs", dest="graphs", help="graphs in pickle format", required=False,
                        default=default_graphs_path)
    parser.add_argument("-a", "--pathways", dest="pathways", help="Pathways list", required=False,
                        default=default_pathways_path)
    parser.add_argument("-n", "--names", dest="names", help="Pathway names", required=False,
                        default=default_pathways_names_path)
    parser.add_argument("-c", "--classes", dest="classes", help="Pathway classes", required=False,
                        default=default_pathways_class_path)

    parser.add_argument("-o", "--outname", dest="outname", help="first part of ouput name", default="summary.kegg")
    parser.add_argument("-w", "--include-weights", dest="include_weights", help="add weights for each KO in output",
                        action='store_true')
    parser.add_argument("-p", "--plot-pathways", dest="plot_pathways", help="Create images with pathways completeness",
                        action='store_true')
    return parser


def load_pathways_data(path_to_graphs, path_to_graphs_names, path_to_graphs_classes):
    """
    Function loads dict of graph that was saved by docker container to graphs.pickle
    :param outdir: path to file with graphs
    :return: dict of graphs
             dict of names of pathways
             dict of classes of pathways
    """
    if os.path.exists(path_to_graphs):
        graph_file = open(path_to_graphs, 'rb')
        graphs = pickle.load(graph_file)
    else:
        logging.error(f'No graphs file found in {path_to_graphs}')
        sys.exit(1)

    pathway_names = {}
    if os.path.exists(path_to_graphs_names):
        with open(path_to_graphs_names, 'r') as file_names:
            for line in file_names:
                line = line.strip().split(':')
                pathway_names[line[0]] = line[1]
    else:
        logging.error(f'No pathways_names file found in {path_to_graphs_names}')
        sys.exit(1)

    pathway_classes = {}
    if os.path.exists(path_to_graphs_classes):
        with open(path_to_graphs_classes, 'r') as file_classes:
            for line in file_classes:
                line = line.strip().split(':')
                pathway_classes[line[0]] = line[1]
    else:
        logging.error(f'No pathways_names file found in {path_to_graphs_classes}')
        sys.exit(1)
    logging.info('Graphs data loaded')
    return graphs, pathway_names, pathway_classes


def get_list_items(input_path, input_list):
    """
    Function creates a list of items that were found by HMMScan
    :param input_path: file with contigs and their KEGG annotations
    :return: list of unique KOs, { contig_name1: [KO1, KO2,...], contig_name2: [...], ...}
    """
    if not input_path and not input_list:
        logging.error("No necessary input table provided")
        sys.exit(1)
    items = []
    dict_KO_by_contigs = {}
    if input_path:
        if os.path.exists(input_path):
            with open(input_path, 'r') as file_in:
                for line in file_in:
                    line = line.strip().split('\t')
                    name = line[0]
                    if name not in dict_KO_by_contigs:
                        dict_KO_by_contigs[name] = []
                    dict_KO_by_contigs[name] += line[1:]
                    items += line[1:]
        else:
            logging.error(f"No file {input_path}")
            sys.exit(1)
    elif input_list:
        if os.path.exists(input_list):
            name = os.path.basename(input_list)
            with open(input_list, 'r') as f:
                list_kos = f.read().strip().split(',')
                if len(list_kos) == 0:
                    logging.error(f"No KOs found in {input_list}")
                else:
                    dict_KO_by_contigs[name] = list_kos
                    items = list_kos
        else:
            logging.error(f"No file {input_list}")
            sys.exit(1)
    else:
        logging.error("No KOs provided")
    logging.info('KOs loaded')
    return list(set(items)), dict_KO_by_contigs


def intersection(lst1, lst2):
    """
    Intersection of two sets
    :param lst1: first input list
    :param lst2: second input list
    :return: intersection in list format
    """
    return list(set(lst1) & set(lst2))


def finding_paths(G):
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


def calculate_percentage(graph, dict_edges, unnecessary_nodes, edges):
    """
    Function returns the percentage of matches of set of edges and graph.
    Example:
            Pathway: A B C. Edges: A -> percentage = 33
    :param graph: input graph of pathway
    :param dict_edges: dict of edges in graph by labels
    :param unnecessary_nodes: list of all nodes
    :param edges: set of nodes
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
    paths_nodes, paths_labels, metrics, indexes_min = finding_paths(graph)
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


def sort_out_pathways(graphs, edges, pathway_names, pathway_classes,
                      contig_name, file_out_summary, weights_of_KOs, include_weights):
    """
    Function sorts out all pathways and prints info about pathway that percentage of intersection more than 0
    :param graphs: Dict of graphs
    :param edges: list of items to intersect with pathways
    :param contig_name == name of contig, or '' for full summary
    :return: -
    """
    dict_sort_by_percentage = {}
    for name_pathway in graphs:
        graph = graphs[name_pathway]
        if intersection(graph[1], edges) == []:
            continue
        else:
            percentage, number_paths, matching_labels, missing_labels = \
                calculate_percentage(graph=graph[0], dict_edges=graph[1], unnecessary_nodes=graph[2], edges=edges)
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
            if include_weights:
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
            output_line = '\t'.join([out_name_pathway, str(percentage), pathway_names[name_pathway],
                                    pathway_classes[name_pathway], matching_current, missing_current])
            file_out_summary.write(output_line + '\n')
    """
    file_out_summary.write('\n******* REMINDER ********')
    file_out_summary.write('Number of nodes: ' + str(len(edges)) + '\n')
    file_out_summary.write('Set of nodes: ' + str(edges) + '\n')
    """


def set_headers(file_summary, contig):
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


def get_weights_for_KOs(graphs):
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

def main():
    parser = parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        graphs, pathway_names, pathway_classes = load_pathways_data(args.graphs, args.names, args.classes)
        edges, dict_KO_by_contigs = get_list_items(args.input_file, args.input_list)
        name_output = args.outname + '.summary.kegg'

        # COMMON INFO
        logging.info('Generating completeness for whole list of KOs...')
        using_graphs = copy.deepcopy(graphs)
        name_output_summary = name_output + '_pathways.tsv'
        file_out_summary = open(name_output_summary, "wt")
        set_headers(file_out_summary, False)
        weights_of_KOs = get_weights_for_KOs(using_graphs)
        sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, '', file_out_summary, weights_of_KOs, args.include_weights)
        file_out_summary.close()
        logging.info('...Done')

        if args.plot_pathways:
            logging.info('Plot pathways images')
            pathways_schema = {}
            if os.path.exists(args.pathways):
                with open(args.pathways, 'r') as pathways_file:
                    for line in pathways_file:
                        fields = line.strip().split(':')
                        pathways_schema[fields[0]] = fields[1]
            else:
                logging.error('No pathways file found')
            pathways = parse_input(input_file=name_output_summary)
            print(pathways)
            plot_graphs(pathways, graphs, pathways_schema)
            logging.info('...Done. Results are in pathways_plots folder')

        # BY CONTIGS
        logging.info('Generating completeness for contigs...')
        name_output_summary = name_output + '_contigs.tsv'
        file_out_summary = open(name_output_summary, "wt")
        set_headers(file_out_summary, True)
        for contig in dict_KO_by_contigs:
            using_graphs = copy.deepcopy(graphs)
            edges = dict_KO_by_contigs[contig]
            sort_out_pathways(using_graphs, edges, pathway_names, pathway_classes, contig, file_out_summary, weights_of_KOs, args.include_weights)
        file_out_summary.close()
        logging.info('...Done')
        logging.info('Bye!')

if __name__ == "__main__":
    main()
