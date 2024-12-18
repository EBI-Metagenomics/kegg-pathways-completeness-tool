## Generate required KEGG modules files

Those instructions contain information how to generate `all_pathways.txt`, `all_pathways_class.txt` and `all_pathways_names.txt`

### Use python script

```commandline
python3 kegg_pathways_completeness/bin/fetch_modules_list.py 
```
Arguments: \
`-o` - output directory name \
`--check-new-modules` - use that flag if you want to check what new modules are missing in current `pathways_data/all_pathways.txt` file in comparison with KEGG API list. Script fetches a list of modules and compares it with currently using modules file. It produces file `new_modules.txt` with list of missing modules.

### Alternative way: bash commands
```bash
mkdir pathways
cd pathways

# copy  http://rest.kegg.jp/list/module 
wget -O - http://rest.kegg.jp/list/module > list_modules.txt

# create names file
cat list_modules.txt | tr '\t' ':'  > all_pathways_names.txt

cat list_modules.txt | cut -f1 > list_pathways.txt

# get DEFINITION file with pathway from KEGG file 
cat list_pathways.txt | while read line; do wget -O - http://rest.kegg.jp/get/$line | grep ^DEFINITION | cut -c 13-; done  > all_pathways_kos.txt
paste list_pathways.txt all_pathways_kos.txt  | tr '\t' ':' > all_pathways.txt

# take CLASS
cat list_pathways.txt | while read line; do wget -O - http://rest.kegg.jp/get/$line | grep ^CLASS | cut -c 13-; done  > all_pathways_class_kos.txt
paste list_pathways.txt all_pathways_class_kos.txt  | tr '\t' ':' > all_pathways_class.txt

# rm tmp files
rm list_modules.txt list_pathways.txt all_pathways_kos.txt all_pathways_class_kos.txt
```

### Notes:
- make sure that each pathway has only KO-s but not other MO-s

# How to generate graphs and plots

1. Generate `graphs.pkl`
```commandline
python3 kegg_pathways_completeness/bin/make_graphs.py \
  -i kegg_pathways_completeness/pathways_data/all_pathways.txt \
  -o kegg_pathways_completeness/pathways_data
```

2. Generate plots. 

Script will generate each module/graph in `.dot` format and save into **dots** folder. Then it will generate `.png` for each `.dot` file and save into **png** folder.
```commandline
python3 kegg_pathways_completeness/bin/make_graphs/generate_schematic_plots.py \
   -l kegg_pathways_completeness/pathways_data/all_pathways.txt \
   -g kegg_pathways_completeness/pathways_data/graphs.pkl \
   -o kegg_pathways_completeness/pathways_data/plots
```

