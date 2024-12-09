# How to generate graphs 

1. Generate `graphs.pkl`
```commandline
python3 kegg_pathways_completeness/bin/make_graphs/make_graphs.py -i kegg_pathways_completeness/pathways_data/all_pathways.txt

mv graphs.pkl graphs/graphs.pkl
```

2. Generate dot-files
```commandline
# generate dot files
python3 kegg_pathways_completeness/bin/make_graphs/get_dot.py -g kegg_pathways_completeness/graphs/graphs.pkl -l kegg_pathways_completeness/pathways_data/all_pathways.txt
```

3. Generate plots for each graph
```commandline
python3 kegg_pathways_completeness/bin/make_graphs/plot.py -l kegg_pathways_completeness/pathways_data/all_pathways.txt
```

4. Put files into correct folder
```commandline
mv dots png kegg_pathways_completeness/graphs/
```

