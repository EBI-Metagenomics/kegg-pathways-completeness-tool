import os
import numpy as np
import argparse
import sys
import pdb
import pickle
import networkx as nx


def download_pathways(outdir):
    """
    Function loads dict of graph that was saved by docker container to graphs.pickle
    :param outdir: path to file with graphs
    :return: dict of graphs
             dict of names of pathways
             dict of classes of pathways
    """
    path_to_graphs = os.path.join(outdir, "graphs.pkl")
    graph_file = open(path_to_graphs, 'rb')
    graphs = pickle.load(graph_file)

    pathway_names = {}
    path_to_graphs_names = os.path.join(outdir, "pathways/all_pathways_names.txt")
    with open(path_to_graphs_names, 'r') as file_names:
        for line in file_names:
            line = line.strip().split(':')
            pathway_names[line[0]] = line[1]

    pathway_classes = {}
    path_to_graphs_classes = os.path.join(outdir, "pathways/all_pathways_class.txt")
    with open(path_to_graphs_classes, 'r') as file_classes:
        for line in file_classes:
            line = line.strip().split(':')
            pathway_classes[line[0]] = line[1]

    return graphs, pathway_names, pathway_classes


def get_list_items(input_path):
    """
    Function creates a list of items that were found by HMMScan
    :param input_path: file with contigs and their KEGG annotations
    :return: list of items
    """
    items = []
    with open(input_path, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            items += line[1:]
    return list(np.unique(items))


def intersection(lst1, lst2):
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
        #print(metrics[len(metrics)-1], paths_nodes[num], paths_labels[num])
    indexes_min = [index for index in range(len(metrics)) if metrics[index] == min(metrics)]
    return paths_nodes, paths_labels, weights, new_weights, indexes_min


def calculate_percentage(graph, dict_edges, edges, name_pathway):
    """
    Function returns the percentage of matches of set of edges and graph.
    Example:
            Pathway: A B C. Edges: A -> percentage = 33
    :param graph: input graph of pathway
    :param dict_edges: dict of edges in graph by labels
    :param edges: set of nodes
    :return: percentage [0:100]
    """
    #print('**********************************************')
    # set weights
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
                #if graph[start][finish][index]['weight'] == 0:  # UNnecessary node
                #    graph[start][finish][index]['weight'] = 1
                graph[start][finish][index]['weight_new'] = 0

    # find the best path(s)
    paths_nodes, paths_labels, weights, new_weights, indexes_min = finding_paths(graph)
    #print('**********************************************')


    percentage = round((1 - 1. * new_weights[num] / weights[num]) * 100, 2)
    if percentage > 0:
        matching_set = set()
        missing_set = set()
        for num in indexes_min:
            #print('==> Path' + str(paths_labels[num]))
            new_labels = paths_labels[num]

            missing_labels = set(new_labels).difference(set(edges))
            missing_set = missing_set.union(missing_labels)
            #print('Missing labels: ', list(missing_labels))

            existing_labels = set(new_labels).intersection(set(edges))
            matching_set = matching_set.union(existing_labels)
            # print('Existing labels: ', set(new_labels).intersection(set(edges)))

            #extra_labels = set(edges).difference(set(new_labels))
            #print('Extra labels: ', extra_labels)
        return percentage, len(indexes_min), ', '.join(list(matching_set)), ', '.join(list(missing_set))

    else:
        return [None for _ in range(4)]
        #print('PATHWAY: ' + name_pathway)
        #print('Found ' + str(len(indexes_min)) + ' paths in PATHWAY ' + name_pathway)
        #print('Percentage = ' + str(percentage))





def sort_out_pathways(graphs, edges, pathway_names, pathway_classes, outdir):
    """
    Function sorts out all pathways and prints info about pathway that percentage of intersection more than 0
    :param graphs: Dict of graphs
    :param edges: list of items to intersect with pathways
    :return: -
    """
    name_output_summary = os.path.join('summary_pathways.txt')
    name_output_matching = os.path.join('matching_ko_pathways.txt')
    name_output_missing = os.path.join('missing_ko_pathways.txt')

    dict_sort_by_percentage = {}
    for name_pathway in graphs:
        graph = graphs[name_pathway]
        if intersection(graph[1], edges) == []:
            continue
        else:
            percentage, kol_paths, matching_labels, missing_labels = \
                calculate_percentage(graph[0], graph[1], edges, name_pathway)
            if percentage != None:
                if percentage not in dict_sort_by_percentage:
                    dict_sort_by_percentage[percentage] = {}
                dict_sort_by_percentage[percentage][name_pathway] = [kol_paths, matching_labels, missing_labels]
    # output Summary
    with open(name_output_summary, 'w') as file_out_summary:
        for percentage in sorted(list(dict_sort_by_percentage.keys()), reverse=True):
            file_out_summary.write('**********************************************\n')
            file_out_summary.write('Percentage = ' + str(percentage) + '\n')
            for name_pathway in dict_sort_by_percentage[percentage]:
                file_out_summary.write(pathway_classes[name_pathway] + ': ' + pathway_names[name_pathway] +
                                       ', ' + name_pathway + '\n')

        file_out_summary.write('\n******* REMINDER ********')
        file_out_summary.write('Number of nodes: ' + str(len(edges)) + '\n')
        file_out_summary.write('Set of nodes: ' + str(edges) + '\n')

    # output matching KOs
    with open(name_output_matching, 'w') as file_out_matching:
        for percentage in sorted(list(dict_sort_by_percentage.keys()), reverse=True):
            file_out_matching.write('**********************************************\n')
            file_out_matching.write('Percentage = ' + str(percentage) + '\n')
            for name_pathway in dict_sort_by_percentage[percentage]:
                file_out_matching.write('\nPathway: ' + name_pathway + '\n')
                file_out_matching.write(dict_sort_by_percentage[percentage][name_pathway][1] + '\n')

    # output missing KOs
    with open(name_output_missing, 'w') as file_out_missing:
        for percentage in sorted(list(dict_sort_by_percentage.keys()), reverse=True):
            file_out_missing.write('**********************************************\n')
            file_out_missing.write('Percentage = ' + str(percentage) + '\n')
            for name_pathway in dict_sort_by_percentage[percentage]:
                file_out_missing.write('\nPathway: ' + name_pathway + '\n')
                missing_str = dict_sort_by_percentage[percentage][name_pathway][2]
                if missing_str == "":
                    missing_str = 'no missing labels for this pathway'
                file_out_missing.write(missing_str + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates Graphs for each contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Each line = pathway", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir",
                        help="Relative path to directory where you want the output file to be stored (default: cwd)",
                        default=".")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        graphs, pathway_names, pathway_classes = download_pathways(args.outdir)
        edges = get_list_items(args.input_file)
        sort_out_pathways(graphs, edges, pathway_names, pathway_classes, args.outdir)
