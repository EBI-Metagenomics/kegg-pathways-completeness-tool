import os
import numpy as np
import argparse
import sys
import pdb
import pickle
import networkx as nx

def set_order_separators(dict_levels):
    """
    Order of parsing: [',', ' ', '+', '-']
    Function returns the parsing separators in the right order according to levels
    ??? not optimal ???
    """
    keys = sorted(list(dict_levels.keys()))
    if keys != []:
        min_level, max_level = [int(keys[0].split('_')[0]), int(keys[len(keys)-1].split('_')[0])+1]
        orders = [str(j)+'_'+i for j in range(min_level, max_level) for i in [',', ' ', '+', '-']]
        new_order = [element for element in orders if element in keys]
    else:
        new_order = []
    return new_order


def add_to_dict_of_levels(dict_levels, c, cur_level, index):
    """
    Function returns the dict of positions according to the level of space or comma
    Example: {'1_,': [14], '2_,': [9], '0_ ': [3], '1_ ': [12]}
        comma of level 1: position 14
        comma of level 2: position 9
        space of level 0: position 3
        space of level 1: position 12
    """
    symbol = str(cur_level) + '_' + c
    if symbol not in dict_levels:
        dict_levels[symbol] = []
    dict_levels[symbol].append(index)
    return dict_levels


def set_brackets(pathway):
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
        if c == '(':
            num += 1
            cur_open.append(num)
            levels_brackets.append(num)
        elif c == ')':
            levels_brackets.append(cur_open[len(cur_open)-1])
            cur_open.pop()
        else:
            levels_brackets.append(-1)
    return levels_brackets


def set_levels(pathway):
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
        if c == ' ' or c == ',' or c == '-' or c == '+':
            dict_levels = add_to_dict_of_levels(dict_levels, c, cur_level, index)
        elif c == '(':
            cur_level += 1
        elif c == ')':
            cur_level -= 1
        else:
            index += 1
            if index < L:
                while pathway[index] not in [' ', ',', '(', ')', '-', '+']:
                    index += 1
                    if index >= L: break
                index -= 1
        index += 1

    return dict_levels


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def check_brackets(pathway, levels_brackets):
    """
    Function checks is this expression in brackets. Returns without if true
    Example: input (A B C)
            return: A B C
    :param pathway: input string expression
    :return: output string expression
    """
    L = len(pathway)
    if pathway[0] == '(' and pathway[L-1] == ')' and levels_brackets[0] == levels_brackets[L-1]:  # check brackets
        return pathway[1:L-1]
    else:
        return pathway


def recursive_parsing(G, dict_edges, unnecessary_nodes, expression, start_node, end_node, weight):
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
    #print(expression, start_node, end_node)

    if expression == '--':  # case --
        name_missing = 'K00000'
        #print('MAKE EDGE: ' + name_missing, start_node, end_node)
        G.add_edge(start_node, end_node, label=name_missing, weight=0, weight_new=0, name='-')
        unnecessary_nodes.append(name_missing)
        if name_missing not in dict_edges:
            dict_edges[name_missing] = []
        dict_edges[name_missing].append([start_node, end_node])
        return G, dict_edges, unnecessary_nodes

    expression = check_brackets(expression, set_brackets(expression))  # delete brackets (expression)
    cur_dict_levels = set_levels(expression)  # define levels of each part of expression
    separators_order = set_order_separators(cur_dict_levels)  # set orders of existing separators
    cur_weight = weight

    if len(separators_order) == 1:  # case: -K....
        if separators_order[0] == '0_-' and expression[0] == '-':
            #print('MAKE EDGE: ' + expression[1:], start_node, end_node)
            G.add_edge(start_node, end_node, label=expression[1:], weight=0, weight_new=0, name='-')
            unnecessary_nodes.append(expression[1:])
            if expression[1:] not in dict_edges:
                dict_edges[expression[1:]] = []
            dict_edges[expression[1:]].append([start_node, end_node])
            return G, dict_edges, unnecessary_nodes

    if separators_order != []:
        # separator
        field = separators_order[0]
        symbol = field.split('_')[1]

        if symbol == '+' or symbol == ' ':
            cur_weight = cur_weight/(len(cur_dict_levels[field])+1)

        separators = list(np.array(sorted(cur_dict_levels[field])))
        cur_sep = 0
        cur_start_node = start_node
        cur_end_node = end_node

        for separator, num in zip(separators, range(len(separators))):

            if symbol == ' ' or symbol == '+' or symbol == '-':
                cur_end_node = len(list(G.nodes()))
                G.add_node(cur_end_node)
            if symbol == '-' and num > 0:
                cur_weight = 0

            subexpression = expression[cur_sep:separator]

            G, dict_edges, unnecessary_nodes = recursive_parsing(G=G,
                                              dict_edges=dict_edges,
                                              unnecessary_nodes=unnecessary_nodes,
                                              expression=subexpression,
                                              start_node=cur_start_node,
                                              end_node=cur_end_node,
                                              weight=cur_weight)
            cur_sep = separator + 1
            if symbol == ' ' or symbol == '+' or symbol == '-':
                cur_start_node = cur_end_node

        num += 1
        if symbol == ' ' or symbol == '+' or symbol == '-':  # nodes and edges
            cur_start_node = cur_end_node
            cur_end_node = end_node
        if symbol == '-' and num > 0:  # weight
            cur_weight = 0

        G, dict_edges, unnecessary_nodes = recursive_parsing(G=G,
                                          dict_edges=dict_edges,
                                          unnecessary_nodes=unnecessary_nodes,
                                          expression=expression[cur_sep:len(expression)],
                                          start_node=cur_start_node,
                                          end_node=cur_end_node,
                                          weight=cur_weight)
        return G, dict_edges, unnecessary_nodes
    else:
        #print('MAKE EDGE: ' + expression, start_node, end_node)
        if cur_weight == 0:
            G.add_edge(start_node, end_node, label=expression, weight=cur_weight, weight_new=cur_weight, name='-')
            unnecessary_nodes.append(expression)
        else:
            G.add_edge(start_node, end_node, label=expression, weight=cur_weight, weight_new=cur_weight, name='node')
        if expression not in dict_edges:
            dict_edges[expression] = []
        dict_edges[expression].append([start_node, end_node])
        return G, dict_edges, unnecessary_nodes


def pathways_processing(input_file, outdir):
    """
    Main function for processing each pathway.
    All pathways were written in one file by lines in format: <name>:<pathway>.
    Function creates dictionary key: name; value: (graph, dict_edges)
    :param input_file: input file with pathways
    :return:
    """
    graphs = {}  # dict of all graphs and their dict of edges (in tuple format)
    with open(input_file, 'r') as file_in:
        for line in file_in:
            line = line.strip().split(':')
            pathway = line[1]
            name = line[0]
            print(name)
            # Graph creation:
            Graph = nx.MultiDiGraph()
            Graph.add_node(0, color='green')
            Graph.add_node(1, color='red')
            # Parsing
            Graph, dict_edges, unnecessary_nodes = recursive_parsing(
                                                  G=Graph,
                                                  dict_edges={},
                                                  unnecessary_nodes=[],
                                                  expression=pathway,
                                                  start_node=0, end_node=1,
                                                  weight=1)
            # Saving
            graphs[name] = tuple([Graph, dict_edges, unnecessary_nodes])
        print('done')
    path_output = os.path.join(outdir, "graphs.pkl")
    f = open(path_output, "wb")
    pickle.dump(graphs, f)
    f.close()


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
        pathways_processing(args.input_file, args.outdir)
