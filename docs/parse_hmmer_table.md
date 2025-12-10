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

You can use hmmscan or hmmsearch

### - hmmscan

```commandline
hmmscan --noali --cut_ga --domtblout <OUTPUT FILE NAME> ${REFDB}/db_kofam/db_kofam.hmm <PROTEIN_SEQS>
 
# parse table
python3 make_hmmer_table_tab_separated.py -i <OUTPUT FILE NAME> -o <output_name for tab-separated file>

python3 parse_hmmer_tab_separated_table.py -t hmmscan -i <output_name for tab-separated file> -f <PROTEIN_SEQS>
```

### - hmmsearch
```commandline
hmmsearch --noali --cut_ga --domtblout <OUTPUT FILE NAME> ${REFDB}/db_kofam/db_kofam.hmm <PROTEIN_SEQS>
 
# parse table
python3 make_hmmer_table_tab_separated.py -i <OUTPUT FILE NAME> -o <output_name for tab-separated file>

python3 parse_hmmer_tab_separated_table.py -t hmmsearch -i <output_name for tab-separated file> -f <PROTEIN_SEQS>
```


