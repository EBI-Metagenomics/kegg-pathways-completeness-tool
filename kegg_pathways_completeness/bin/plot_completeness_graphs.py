#!/usr/bin/env python3

import argparse
import sys
import pickle
import networkx as nx
import logging
import os
import graphviz

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

def parse_input(input_file):
    pathways = {}
    if os.path.exists(input_file):
        with open(input_file, 'r') as file_in:
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
    else:
        logging.error(f'File {input_file} does not exist')
    return pathways


def plot_graphs(pathways, graphs, pathways_schema):
    for name in pathways:
        graph = graphs[name]
        logging.info(f'Plotting {name}')
        presented_ko = pathways[name].split(',')
        dot = create_dot(name, presented=presented_ko, graph=graph, pathways_schema=pathways_schema[name])
        dot.render(directory='pathways_plots', filename=name, format='png')
        #Image(f'kegg/{name}.png')

def create_dot(name, presented, graph, pathways_schema):
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
            if 1/weight > max_weight:
                max_weight = int(1/weight)
        if weight == 1 or weight == 0:
            weight_str = str(weight)
        else:
            weight_str = '1/' + str(int(1/weight))
        color = 'red' if label in presented else 'black'
        dot.edge(str(from_node), str(to_node), label=label + ' \n [' + weight_str + ']', color=color)
    return dot

def main():
    parser = argparse.ArgumentParser(description="Script generates Graphs for each contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Output with completeness", required=True)
    parser.add_argument("-g", "--graphs", dest="graphs", help="graphs in pickle format", required=False,
                        default="graphs/graphs.pkl")
    parser.add_argument("-p", "--pathways", dest="pathways", help="Pathways of kos", required=False,
                        default="pathways_data/all_pathways.txt")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        pathways_schema = {}
        if os.path.exists(args.pathways):
            with open(args.pathways, 'r') as pathways_file:
                for line in pathways_file:
                    fields = line.strip().split(':')
                    pathways_schema[fields[0]] = fields[1]
        else:
            logging.error('No pathways file found')

        if os.path.exists(args.graphs):
            graphs = pickle.load(open(args.graphs, 'rb'))
            plot_graphs(parse_input(input_file=args.input_file), graphs, pathways_schema)
        else:
            logging.error('No graphs file found')

if __name__ == "__main__":
    main()

