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