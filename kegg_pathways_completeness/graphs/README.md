# How to generate graphs 

1. Generate `graphs.pkl`
```commandline
python3 bin/make_graphs/make_graphs.py -i pathways_data/all_pathways.txt

mv graphs.pkl graphs/graphs.pkl
```

2. Generate dot-files
```commandline
# generate dot files
python3 bin/make_graphs/get_dot.py -g graphs/graphs.pkl -l pathways_data/all_pathways.txt
```

3. Generate plots for each graph
```commandline
python3 bin/make_graphs/plot.py -l pathways_data/all_pathways.txt
```

4. Put files into correct folder
```commandline
mv dots png graphs/
```

