# kegg-pathways-completeness tool

The tool counts completeness of each KEGG pathway for protein sequence. 

Please read **Theory** section with detailed explanation in the bottom of README. 

Current version of pathways saved into **[pathways_data](pathways_data)** and graphs were [pre-built](graphs/README.md) and [saved](graphs/updates/pipeline-v5/graphs.pkl) into pkl format. 

`*Pipeline-v5 data has 394 modules.*`

`Updated (from 07/03/2024) data in this repo has 481 modules.`

Previous updates:
- 27/04/2023 has 475 modules.

If you need to update existing pathways data and graphs follow this [instruction](pathways_data/README.md).

These files are also available on EBI MGnify FTP and can be downloaded using [download.sh](download.sh)

## Calculate pathways completeness
This script requires [hmmsearch table](tests/test_data/test-input/test) run on KEGG profiles with annotated sequences (preferable) **OR** [file with list](tests/test_data/test-input/test_list.txt) of KOs.
If you don't have this table follow [instructions](src/README.md) how to generate it first.

#### Run using conda 
```commandline
conda create --name kegg-env
conda activate kegg-env

pip3 install requirements.txt

export INPUT='tests/test_data/test-input/test'  # path to hmm-result table
export OUTPUT='test-out'  # prefix for output

# hmmtable as input
python3 bin/give_pathways.py \
  -i ${INPUT} \
  -o ${OUTPUT}

# KOs list as input
python3 bin/give_pathways.py \
  -l 'tests/test_data/test-input/test_list.txt' \
  -o ${OUTPUT}
```
Check example of output [here](tests/test_data/test-output). \
`kegg_pathways.tsv` has pathways completeness calculated by all KOs in given input file \
`kegg_contigs.tsv` has pathways completeness calculated per each contig (first column contains name of contig).


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
**NOTE**: please install graphviz \
If you want to see what edges were chosen to complete the graph of completeness you can plot them adding **_--plot-pathways_** argument. \
```commandline
python3 bin/give_pathways.py -i ${INPUT} -o ${OUTPUT} --plot-pathways
```
You can also run plotting script separately:
```commandline
python3 bin/plot_completeness_graphs.py -i output_with_pathways_completeness
```

Example,

![M00050.png](tests%2Ftest_data%2Ftest-output%2Fplots%2FM00050.png)

more examples for test data [here](tests/test_data/test-output/plots)


## Theory: 
### Pathways to graphs 
KEGG provides a representation of each pathway as specific expression of KOs.
ex: **A ((B,C) D,E) (A+F)** \
where A, B, C, D, E, F are KOs \
**space** means AND \
**comma** means OR \
**plus** means essential component \
**minus** means optional component
Each expression was recursively [converted](bin/make_graphs/make_graphs.py) into directed graph using NetworkX. First node has number 0 and the last number 1. Each edge corresponds to KO. 

![ex1.png](src%2Fimg%2Fex1.png)

### Completeness
In order to count pathways completeness each graph was made weighted. Default weight of each edge is 0. \
Let's imagine there is a set of KOs predicted by annotation. If KO is presented in pathway - corresponding edge receives weight = 1 (or 0 if edge is optional or another value if edge is connected by +). \
After that [script](bin/give_pathways.py) searches the most weighted path from node 0 to node 1 (`graph_weight`). 
`max_graph_weight` calculated in assumption all KOs are presented. \
``
completeness = graph_weight/max_graph_weight * 100%
``

![ex2.png](src%2Fimg%2Fex2.png)


## Create plots for all pathways
There are [plots](graphs/png) for every pathway as graph representation.
If you need to re-generate them follow [instruction](graphs/README.md).

