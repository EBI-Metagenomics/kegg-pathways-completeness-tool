# kegg-pathways-completeness tool

The tool counts completeness of each KEGG pathway for protein sequence. 
Current version of pathways saved into **[pathways_data](pathways_data)** and graphs were [pre-built](graphs/README.md) and [saved](graphs/graphs.pkl) into pkl format. 
These files are available on EBI MGnify FTP and can be downloaded using [download.sh](download.sh)

## Calculate pathways completeness
This script requires hmmsearch table run on KEGG profiles with annotated sequences.
If you don't have this table follow [instructions](src/README.md) how to generate it.

#### Run using env 
```commandline
export INPUT="path to hmm-result table"
export OUTPUT=result
pip3 install requirements.txt
python3 Tools/give_pathways.py \
  -i ${INPUT} \
  -g help_files/graphs.pkl \
  -c help_files/all_pathways_class.txt \
  -n help_files/all_pathways_names.txt \
  -o ${OUTPUT}
```

#### Run using docker
Results can be fould be in folder "results". Final annotated pathways would be in folder "results/pathways
```commandline
export INPUT="path to hmm-result table"
docker \
    run \
    -i \
    --workdir=/results \
    --volume=`pwd`/results:/results:rw \
    --volume=${INPUT}:/files/input_table.tsv:ro \
    kegg_pathways:latest \
    /tools/run_pathways.sh \
    -i /files/input_table.tsv
```

## Update existing pathways data
If you need to update existing pathways data and graphs follow this [instruction]().

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


## Create plots of pathways
There are [plots](graphs/png) for every pathway as graph representation.
If you need to re-generate them follow [instruction](graphs/README.md).

