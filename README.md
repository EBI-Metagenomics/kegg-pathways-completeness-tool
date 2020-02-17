# kegg-cwl

## How to run

#### Terminal
```bash
# get help_files
mkdir help_files
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/graphs.pkl.gz -P help_files && gzip -d help_files/graphs.pkl.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_names.txt.gz -P help_files && gzip -d help_files/all_pathways_names.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_class.txt.gz -P help_files && gzip -d help_files/all_pathways_class.txt.gz

#Prepare table in format: accertion \t KO_#1 \t KO_#2 ...
export INPUT=<your table from step 3>
export OUTPUT=pathways_results
 
pip3 install networkx
python3 Tools/give_pathways.py -i ${INPUT} -g help_files/graphs.pkl -c help_files/all_pathways_class.txt -n help_files/all_pathways_names.txt -o ${OUTPUT}
```

#### Docker
```bash
# ! go to folder with DockerFile

export TABLE= ! fill in full path !

# build docker
docker build -t kegg_pathways .

# run docker
docker \
    run \
    -i \
    --workdir=/results \
    --volume=`pwd`/results:/results:rw \
    --volume=${TABLE}:/files/table_file.txt:ro \
    kegg_pathways:latest \
    /tools/run_pathways.sh \
    -i \
    /files/table_file.txt

# Results fould be in folder "results". Final annotated pathways would be in folder "results/pathways":
    file with "_pathways" - results by all contigs
```

### Prepare files from KEGG (all these files could be found in **pathways** folder)
```bash
# download list of all pathways
curl -s http://rest.kegg.jp/list/module | cut -f1 | cut -c 4- > pathways/list_pathways.txt

# get DEFINITION file with pathway from KEGG file (saving in format: <name:pathway>)
cat pathways/list_pathways.txt | while read line; do echo "$line:" && curl -s http://rest.kegg.jp/get/$line | grep ^DEFINITION | cut -c 13-; done |  sed -z 's|:\n|:|g' > pathways/all_pathways.txt
# !!! TODO make sure that each pathway has only KO-s but not other MO-s

parallel -k echo -n '{}:' ';' curl -s http://rest.kegg.jp/get/{} '|' grep ^NAME '|' cut -c 13- :::: pathways/list_pathways.txt > pathways/all_pathways_names.txt
# taking CLASS
parallel -k echo -n '{}:' ';' curl -s http://rest.kegg.jp/get/{} '|' grep ^CLASS '|' cut -c 13- :::: pathways/list_pathways.txt > pathways/all_pathways_class.txt
```

### Create graphs
```bash
python3 graphs/make_graph/make_graphs.py -i pathways/all_pathways.txt
mkdit graphs/dots graphs/png
python3 graphs/get_dot.py -i graphs/make_graph/graphs.pkl
python3 plot.py
```

