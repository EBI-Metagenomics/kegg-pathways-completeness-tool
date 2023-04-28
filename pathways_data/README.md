### Prepare files from KEGG (all these files could be found in **pathways** folder)

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

NOTE: make sure that each pathway has only KO-s but not other MO-s