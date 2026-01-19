# How to Generate KO-Annotated Table

This guide explains how to generate and process HMMER output to extract KEGG Ortholog (KO) annotations for each contig.


### Required Software
- Python 3 with BioPython
- [HMMER](http://hmmer.org/) (hmmscan or hmmsearch)

### Input 
1. **Protein sequences** (FASTA format) - protein sequences from your contigs
2. **KEGG KOfam HMM database** - can be downloaded from [KEGG KOfam](https://www.genome.jp/tools/kofamkoala/) or generated with [MGnify database generation pipeline](https://github.com/EBI-Metagenomics/reference-databases-preprocessing-pipeline) using `--generate_kofam_db`:


### Output 
Tab-separated file with one contig per line:
- Column 1: Contig name (from FASTA file)
- Columns 2+: All KO identifiers assigned to proteins in that contig

Example:
```
contig_1    K00001  K00002  K00003
contig_2    K00010  K00015
contig_3    K00001
```

## Step 1: Run HMMER search

You can use either `hmmscan` or `hmmsearch`. Both will produce similar results, but the column order in the output differs (check [explanation](http://cryptogenomicon.org/hmmscan-vs-hmmsearch-speed-the-numerology.html) for more details).

### Option A: Using hmmscan (recommended)

```bash
hmmscan \
  --noali \
  --cut_ga \
  --domtblout hmmer_output.tbl \
  db_kofam.hmm proteins.faa
```

### Option B: Using hmmsearch

```bash
hmmsearch \
  --noali \
  --cut_ga \
  --domtblout hmmer_output.tbl \
  db_kofam.hmm proteins.faa
```

**Parameters explained:**
- `--noali`: Don't output alignments (saves space and time)
- `--cut_ga`: Use gathering threshold for filtering hits
- `--domtblout`: Output domain hit table format
- `db_kofam.hmm`: KOfam HMM database file
- `proteins.faa`: Your protein sequences in FASTA format

## Step 2: Process HMMER Output

Use the `parse_hmmer_table.py` tool to extract KO annotations per contig. The tool automatically handles the conversion from HMMER's space-separated format to tab-separated format and then parses it to extract KO annotations. The parser supports both input formats: hmmsearch and hmmscan (specify with `--tool` argument). 

### Basic Usage

```bash
parse_hmmer_table \
  -i HMMER_OUTPUT.tbl \
  -f proteins.faa \
  -t tool [hmmscan / hmmsearch] \
  -o ko_annotations.tsv \
  [--save-intermediate hmmer_intermediate.tsv]
```

### Required Parameters
- `-i, --input`: Input HMMER domtblout file (output from Step 1)
- `-f, --fasta`: FASTA file with protein sequences (same file used in HMMER search)
- `-t, --tool`: HMMER tool used, either `hmmscan` or `hmmsearch`
- `-o, --output`: Output file for KO annotations per contig

### Optional Parameters
- `--save-intermediate`: Path to save intermediate tab-separated table (optional)