# How to generate hmm table

```commandline
hmmscan --noali --cut_ga --domtblout <OUTPUT FILE NAME> ${REFDB}/db_kofam/db_kofam.hmm <PROTEINS FILE>
 
# prepare table
python3 hmmscan_tab.py -i <your table> -o <output_name>
python3 parsing_hmmscan.py -i <output_name> -f <your fasta that you used for hmmscan>
```
