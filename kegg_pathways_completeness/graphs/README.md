# How to generate graphs and plots

1. Generate `graphs.pkl`
```commandline
python3 kegg_pathways_completeness/bin/make_graphs/make_graphs.py \
  -i kegg_pathways_completeness/pathways_data/all_pathways.txt \
  -o kegg_pathways_completeness/graphs
```

2. Generate plots. 

Script will generate each module/graph in `.dot` format and save into **dots** folder. Then it will generate `.png` for each `.dot` file and save into **png** folder.
```commandline
python3 kegg_pathways_completeness/bin/make_graphs/generate_schematic_plots.py \
   -l kegg_pathways_completeness/pathways_data/all_pathways.txt \
   -g kegg_pathways_completeness/graphs/graphs.pkl \
   -o kegg_pathways_completeness/graphs
```


