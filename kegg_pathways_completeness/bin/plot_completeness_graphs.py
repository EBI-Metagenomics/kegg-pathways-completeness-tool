#!/usr/bin/env python3

import argparse
import sys
import pickle
import networkx as nx
import logging
import os
import graphviz

logging.basicConfig(encoding='utf-8', level=logging.DEBUG)

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Script generates plots for each contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Output with completeness", required=True)
    parser.add_argument("-g", "--graphs", dest="graphs", help="graphs in pickle format", required=False,
                        default="graphs/graphs.pkl")
    parser.add_argument("-p", "--pathways", dest="pathways", help="Pathways of kos", required=False,
                        default="pathways_data/all_pathways.txt")
    parser.add_argument("-o", "--outdir", dest="outdir", help="Path to output directory", required=False,
                        default="pathways_plots")
    return parser.parse_args(argv)


class PlotModuleCompletenessGraph():
    def __init__(
            self,
            completeness_file:str,
            graphs_pkl: str,
            modules_list: str,
            outdir: str
    ):
        """
        Script generates a plots with colored edges depending on presense/absense.
        If edge is presented - it is colored in red, otherwise - black
        :param completeness_file: output of give_pathways.py
        :param graphs_pkl: graphs.pkl
        :param modules_list: all_pathways.txt
        :param outdir: name of output directory, default=
        """
        self.completeness_file = completeness_file
        self.graphs_pkl = graphs_pkl
        self.modules_list = modules_list
        self.outdir = outdir

    def parse_input(self):
        pathways = {}
        if os.path.exists(self.completeness_file):
            with open(self.completeness_file, 'r') as file_in:
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
            logging.error(f'File {self.completeness_file} does not exist')
        return pathways

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

    def plot_graphs(self, pathways, graphs, pathways_schema):
        for name in pathways:
            graph = graphs[name]
            logging.info(f'Plotting {name}')
            presented_ko = pathways[name].split(',')
            dot = self.create_graph(name, presented=presented_ko, graph=graph, pathways_schema=pathways_schema[name])
            dot.render(directory=self.outdir, filename=name, format='png')
            #Image(f'kegg/{name}.png')

    def generate_plot_for_completeness(self):
        pathways_schema = {}
        if os.path.exists(self.modules_list):
            with open(self.modules_list, 'r') as pathways_file:
                for line in pathways_file:
                    fields = line.strip().split(':')
                    pathways_schema[fields[0]] = fields[1]
        else:
            logging.error('No pathways file found')

        if os.path.exists(self.graphs_pkl):
            with open(self.graphs_pkl, 'rb') as file_graph:
                graphs = pickle.load(file_graph)
                self.plot_graphs(self.parse_input(), graphs, pathways_schema)
        else:
            logging.error('No graphs file found')


def main():
    args = parse_args(sys.argv[1:])
    plot_completeness_generator = PlotModuleCompletenessGraph(
        completeness_file=args.input_file,
        graphs_pkl=args.graphs,
        modules_list=args.pathways,
        outdir=args.outdir,
    )
    plot_completeness_generator.generate_plot_for_completeness()


if __name__ == "__main__":
    main()

