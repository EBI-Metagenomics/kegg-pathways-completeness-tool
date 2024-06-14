# How to generate graphs 

```commandline
# generate graphs.pkl
python3 bin/make_graphs/make_graphs.py -i pathways_data/all_pathways.txt

mv graphs.pkl graphs/graphs.pkl

# generate dot files
python3 bin/make_graphs/get_dot.py -g graphs/graphs.pkl -l pathways_data/all_pathways.txt

python3 bin/make_graphs/plot.py -l pathways_data/all_pathways.txt

mv dots png graphs/
```