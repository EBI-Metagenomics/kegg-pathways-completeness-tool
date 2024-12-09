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

import os
import argparse
import sys
import pickle
import networkx as nx
import subprocess


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates graphs in dot format for each Graph")
    parser.add_argument("-g", "--input", dest="graphs", help="graphs.pkl", required=True)
    parser.add_argument("-l", "--pathways", dest="list_pathways", help="all_pathways.txt", required=True)
    parser.add_argument("-o", "--outdir", dest="outdir", help="Path to output directory", required=False,
                        default=os.getcwd())
    return parser.parse_args(argv)


class PlotGenerator:
    def __init__(
            self,
            input_graphs: str,
            input_list_modules: str,
            outdir: str,
    ):
        """
        Helper to generate graphs in dot format for further parsing into png representation.
        :param input_graphs: graphs.pkl
        :param input_list_modules: all_pathways.txt in format MXXX:KOXX...
        """
        self.input_graphs = input_graphs
        self.input_list_modules = input_list_modules
        self.outdir = outdir
        dot_dir_name = 'dots'
        png_dir_name = 'png'
        self.dots_dir = os.path.join(outdir, dot_dir_name)
        self.plots_dir = os.path.join(outdir, png_dir_name)

    def create_dot(self, graph, name, pathway):
        if not os.path.exists(self.dots_dir):
            os.mkdir(self.dots_dir)
        with open(os.path.join(self.dots_dir, name + '.dot'), 'w') as dot_file:
            dot_file.write("digraph G {\n"
                           "graph [label=\"" + name + "\n" + pathway + "\",fontsize=20];\n"
                           "node [shape=box,style=filled];\n"
                           "edge [len=3,color=grey];\n"
                           "{node [width=.3,height=.3,shape=octagon,style=filled,color=skyblue] ")
            edges = graph[0].edges
            for node in graph[0].nodes:
                dot_file.write(str(node) + ' ')
            dot_file.write('}\n')
            for edge, count in zip(edges, range(len(edges))):
                from_node = edge[0]
                to_node = edge[1]
                number = edge[2]
                label = edges._adjdict[from_node][to_node][number]['label']
                weight = edges._adjdict[from_node][to_node][number]['weight']
                if weight == 1 or weight == 0 :
                    weight_str = str(weight)
                else:
                    weight_str = '1/' + str(int(1/weight))
                dot_file.write(str(from_node) + ' -> ' + str(to_node) + ' [label="' + label + ' [' + weight_str + ']"];\n')
            dot_file.write('}')

    def generate_dot_files(self):
        with open(self.input_graphs, 'rb') as file_graph, open(self.input_list_modules, 'r') as list_pathways:
            graphs = pickle.load(file_graph)
            for line in list_pathways:
                line = line.strip()
                name, pathway = line.split(':')
                print(f"Dot for {name}")
                self.create_dot(graphs[name], name, pathway)

    def generate_png(self, name):
        bashCommand = f"neato -T png {self.dots_dir}/{name}.dot > {self.plots_dir}/{name}.png"
        subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)

    def generate_plots(self):
        if not os.path.exists(self.plots_dir):
            os.mkdir(self.plots_dir)
        with open(self.input_list_modules, 'r') as file_pathways:
            for line in file_pathways:
                name = line.strip().split(':')[0]
                print(f"Plot for {name}")
                self.generate_png(name=name)

def main():
    args = parse_args(sys.argv[1:])
    plot_generator = PlotGenerator(
        input_graphs=args.graphs,
        input_list_modules=args.list_pathways,
        outdir=args.outdir
    )
    plot_generator.generate_dot_files()
    plot_generator.generate_plots()


if __name__ == "__main__":
    main()

