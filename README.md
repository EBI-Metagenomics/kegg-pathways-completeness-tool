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

## Create plots of pathways
There are [plots](graphs/png) for every pathway.
If you need to re-generate them follow [instruction](graphs/README.md).

