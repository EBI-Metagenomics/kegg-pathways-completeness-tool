# How to generate KOs-annotated hmm-table

## Input 
Protein fasta file (PROTEIN_SEQS)

## Required tools
- python3
- [hmmer](http://hmmer.org/)
- reference database: can be downloaded from MGnify FTP:
```
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/db_kofam.hmm*
```

## Run

```commandline
hmmscan --noali --cut_ga --domtblout <OUTPUT FILE NAME> ${REFDB}/db_kofam/db_kofam.hmm <PROTEIN_SEQS>
 
# parse table
python3 make_hmmscan_table_tab_separated.py -i <OUTPUT FILE NAME> -o <output_name for tab-separated file>

python3 parse_hmmscan_tab_separated_table.py -i <output_name for tab-separated file> -f <PROTEIN_SEQS>
```
