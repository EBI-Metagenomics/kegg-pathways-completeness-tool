# kegg-pathways-completeness tool

The tool counts completeness of each KEGG modules pathway for protein sequence. 

Please read **Theory** section with detailed explanation in the bottom of README. 



## Calculate pathways completeness

Current version of tool has 481 KEGG modules updated 07/03/2024. 
If you need to manually update existing pathways and graphs or generate data for not existing pathway - follow this [instruction](kegg_pathways_completeness/pathways_data/README.md).

### Input arguments description

**Required arguments:** 

_interchangeable arguments:_
- input table (`-i`/`--input`): hmmsearch table ([example](tests/fixtures/give_pathways/test_pathway.txt)) that was run on KEGG profiles with annotated sequences (preferable). If you don't have this table follow [instructions](src/README.md) how to generate it.
- file with KOs list (`-l`/`--input-list`): comma separated file with list of KOs ([example](tests/fixtures/give_pathways/test_kos.txt)).

_pathways arguments (all required):_

This repository has a set of required files pre-generated. Current version of data was saved into **[pathways_data](kegg_pathways_completeness/pathways_data)**. 

- list of KEGG modules in KOs notation (`-a`/`--pathways`) (latest [all_pathways.txt](kegg_pathways_completeness%2Fpathways_data%2Fall_pathways.txt))
- list of classes of KEGG modules (`-c`/`--classes`) (latest [all_pathways_class.txt](kegg_pathways_completeness%2Fpathways_data%2Fall_pathways_class.txt))
- list of names of KEGG modules (`-n`/`--names`) (latest [all_pathways_names.txt](kegg_pathways_completeness%2Fpathways_data%2Fall_pathways_names.txt))

_graphs:_

In order to generate graphs all pathways were parsed with networkx library. Every graph is presented in .png format in [png](kegg_pathways_completeness/graphs/png) and .dot format in [dots](kegg_pathways_completeness/graphs/dots). Pathway and weights of each KO can be checked easily with .png image.
Instructions how to build graphs.pkl are [provided](kegg_pathways_completeness/graphs/README.md). 

- graphs constructed from each module (`-g`/`--graphs`) (latest [graphs.pkl](kegg_pathways_completeness%2Fgraphs%2Fgraphs.pkl))

**Optional arguments:**

- output prefix (`-o`/`--outname`): prefix for output tables (`-o test_kos` in [example](tests/fixtures/give_pathways/output/test_kos.summary.kegg_contigs.tsv))
- add weights information to output files (`-w`/`--include-weights`). Output table will have a weight of each KO edge in pathway graph, example K00942(0.25) means that KO has a 0.25 importance in given pathway. Example of [output](tests/fixtures/give_pathways/output/test_weights.summary.kegg_pathways.tsv)
- plot presented pathways (`p`/`--plot-pathways`): PNG contains a schematic representation of pathway. Presented KOs marked with Red edges. Example [M00001](tests/fixtures/give_pathways/output/pathways_plots/M00001.png)

### Example of output

Check example of output [here](tests/fixtures/give_pathways/output). 
- `*kegg_pathways.tsv` has pathways completeness calculated by all KOs in given input file \
- `*kegg_contigs.tsv` has pathways completeness calculated per each contig (first column contains name of contig).
- `pathways_plots` example of plots and graphs generated with `--plot-pathways` argument. 
- `*weights*.tsv` example of output generated with `--include-weights` argument

## Installation
That tool was published in Pypi and Bioconda:
#### using Pip
```commandline
pip install kegg-pathways-completeness
```

## Run
```
give_pathways -l {INPUT_LIST}
```

#### Install from source using venv/conda env
```commandline
conda create --name kegg-env
conda activate kegg-env

pip3 install -r requirements.txt

# Run
# hmmtable as input
python3 kegg_pathways_completeness/bin/give_pathways.py \
  -i 'tests/fixtures/give_pathways/test_pathway.txt' \
  -o test_pathway

# KOs list as input
python3 kegg_pathways_completeness/bin/give_pathways.py \
  -l 'tests/fixtures/give_pathways/test_kos.txt' \
  -o test_list_kos
```

#### Run using docker 
Results can be found in folder `results`. Final annotated pathways would be in folder `results/pathways`
```commandline
export INPUT="path to hmm-result table"
docker \
    run \
    -i \
    --workdir=/results \
    --volume=`pwd`/results:/results:rw \
    --volume=${INPUT}:/files/input_table.tsv:ro \
    quay.io/microbiome-informatics/kegg-completeness:v1.1 \
    /tools/run_pathways.sh \
    -i /files/input_table.tsv
```

## Plot pathways completeness
**NOTE**: please make sure you have [**graphviz**](https://graphviz.org/) installed \
You can also run plotting script separately:
```commandline
python3 kegg_pathways_completeness/bin/plot_completeness_graphs.py -i output_with_pathways_completeness
```

Example,

![M00050.png](tests/fixtures/give_pathways/output/pathways_plots/M00050.png)

more examples for test data [here](tests/fixtures/give_pathways/output/pathways_plots)


## Theory: 
### Pathways to graphs 
KEGG provides a representation of each pathway as specific expression of KOs.
ex: **A ((B,C) D,E) (A+F)** \
where A, B, C, D, E, F are KOs \
**space** means AND \
**comma** means OR \
**plus** means essential component \
**minus** means optional component
Each expression was recursively [converted](kegg_pathways_completeness/bin/make_graphs/make_graphs.py) into directed graph using NetworkX. First node has number 0 and the last number 1. Each edge corresponds to KO. 

![ex1.png](src%2Fimg%2Fex1.png)

### Completeness
In order to count pathways completeness each graph was made weighted. Default weight of each edge is 0. \
Let's imagine there is a set of KOs predicted by annotation. If KO is presented in pathway - corresponding edge receives weight = 1 (or 0 if edge is optional or another value if edge is connected by +). \
After that [script](kegg_pathways_completeness/bin/give_pathways.py) searches the most weighted path from node 0 to node 1 (`graph_weight`). 
`max_graph_weight` calculated in assumption all KOs are presented. \
``
completeness = graph_weight/max_graph_weight * 100%
``

![ex2.png](src%2Fimg%2Fex2.png)

