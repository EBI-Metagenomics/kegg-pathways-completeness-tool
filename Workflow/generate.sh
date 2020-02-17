#!/bin/bash

python3 graphs/make_graph/make_graphs.py -i pathways/all_pathways.txt

python3 graphs/get_dot.py -i graphs/make_graph/graphs.pkl

python3 plot.py