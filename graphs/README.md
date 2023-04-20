# How to generate graphs 

```commandline
# generate graphs.pkl
python3 bin/make_graphs/make_graphs.py -i pathways_data/all_pathways.txt

# generate dot files
python3 bin/make_graphs/get_dot.py -i graphs/make_graph/graphs.pkl

python3 plot.py
```