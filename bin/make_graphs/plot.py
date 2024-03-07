import subprocess
import argparse
import sys
import os

def plot_graph(name):

    #pos = nx.spring_layout(G)
    #write_dot(G, name_graph + '.dot')
    bashCommand = "neato -T png dots/{n}.dot > png/{n}.png".format(n=name)
    subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates Graphs for each contig")
    parser.add_argument("-l", "--pathways", dest="list_pathways", help="all_pathways.txt", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        if not os.path.exists('png'):
            os.mkdir('png')
        with open(args.list_pathways, 'r') as file_pathways:
            for line in file_pathways:
                plot_graph(line.strip().split(':')[0])