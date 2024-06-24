# Updates:

**Previous [updates](kegg_pathways_completeness/graphs/updates):**
- 27/04/2023 has 475 modules.
- MGnify [pipeline-v5](https://github.com/EBI-Metagenomics/pipeline-v5) uses 394 modules.

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