import os
import argparse
import sys
import pickle
import networkx as nx


def create_dot(graph, name, pathway):
    if not os.path.exists('dots'):
        os.mkdir('dots')
    with open(os.path.join('dots', name + '.dot'), 'w') as dot_file:
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates Graphs for each contig")
    parser.add_argument("-g", "--input", dest="input_file", help="graphs.pkl", required=True)
    parser.add_argument("-l", "--pathways", dest="list_pathways", help="all_pathways.txt", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        graph_file = open(args.input_file, 'rb')
        graphs = pickle.load(graph_file)
        with open(args.list_pathways, 'r') as list_pathways:
            for line in list_pathways:
                line = line.strip()
                name, pathway = line.split(':')
                print(name)
                create_dot(graphs[name], name, pathway)
