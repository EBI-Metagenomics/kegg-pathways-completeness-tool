# How to Generate KO-Annotated Table from HMMER Results

This guide explains how to process HMMER output to extract KEGG Ortholog (KO) annotations for each contig.

## Prerequisites

### Required Software
- Python 3 with BioPython
- [HMMER](http://hmmer.org/) (hmmscan or hmmsearch)

### Input Files
1. **Protein sequences** (FASTA format) - protein sequences from your contigs
2. **KEGG KOfam HMM database** - can be downloaded from MGnify FTP:
   ```bash
   wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/db_kofam.hmm.gz
   gunzip db_kofam.hmm.gz
   ```

### Output Format
Tab-separated file with one contig per line:
- Column 1: Contig name (from FASTA file)
- Columns 2+: All KO identifiers assigned to proteins in that contig

Example:
```
contig_1    K00001  K00002  K00003
contig_2    K00010  K00015
contig_3    K00001
```

## Step 1: Run HMMER Search

You can use either `hmmscan` or `hmmsearch`. Both will produce similar results, but the column order in the output differs.

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

Use the `parse_hmmer_table` tool to extract KO annotations per contig. The tool automatically handles the conversion from HMMER's space-separated format to tab-separated format and then parses it to extract KO annotations.

### Basic Usage

```bash
parse_hmmer_table \
  -i hmmer_output.tbl \
  -f proteins.faa \
  -t hmmscan \
  -o ko_annotations.tsv
```

### With Intermediate File (Optional)

If you need to keep the intermediate tab-separated HMMER table for debugging or reuse:

```bash
parse_hmmer_table \
  -i hmmer_output.tbl \
  -f proteins.faa \
  -t hmmscan \
  -o ko_annotations.tsv \
  --save-intermediate hmmer_intermediate.tsv
```

This will save both:
- `hmmer_intermediate.tsv` - intermediate tab-separated HMMER table
- `ko_annotations.tsv` - final KO annotations per contig

## Parameters Explained

### Required Parameters
- `-i, --input`: Input HMMER domtblout file (output from Step 1)
- `-f, --fasta`: FASTA file with protein sequences (same file used in HMMER search)
- `-t, --tool`: HMMER tool used, either `hmmscan` or `hmmsearch`
- `-o, --output`: Output file for KO annotations per contig

### Optional Parameters
- `--save-intermediate`: Path to save intermediate tab-separated table (optional)

**Important:** The `-t` parameter must match the HMMER tool you used in Step 1, as the column order differs between hmmscan and hmmsearch.

## Complete Example Workflow

Here's a complete example using hmmscan:

```bash
# 1. Download KOfam database (if not already available)
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/db_kofam.hmm.gz
gunzip db_kofam.hmm.gz

# 2. Run HMMER search
hmmscan \
  --noali \
  --cut_ga \
  --domtblout my_proteins.hmmer.tbl \
  db_kofam.hmm my_proteins.faa

# 3. Process HMMER output to get KO annotations per contig
parse_hmmer_table \
  -i my_proteins.hmmer.tbl \
  -f my_proteins.faa \
  -t hmmscan \
  -o ko_annotations.tsv

# Output file: ko_annotations.tsv
```

## Example with hmmsearch

```bash
# 1. Run HMMER search with hmmsearch
hmmsearch \
  --noali \
  --cut_ga \
  --domtblout my_proteins.hmmer.tbl \
  db_kofam.hmm my_proteins.faa

# 2. Process with correct tool parameter
parse_hmmer_table \
  -i my_proteins.hmmer.tbl \
  -f my_proteins.faa \
  -t hmmsearch \
  -o ko_annotations.tsv
```

## Example with Output in Subdirectory

The tool automatically creates the output directory if it doesn't exist:

```bash
parse_hmmer_table \
  -i my_proteins.hmmer.tbl \
  -f my_proteins.faa \
  -t hmmscan \
  -o results/ko_annotations.tsv
```

## Troubleshooting

### "Contig not found in FASTA file" warnings
This can happen if contig names in the HMMER output differ from the FASTA file. The script tries to match partial names, but ensure your FASTA headers are consistent.

### Empty output file
- Check that HMMER found hits (review the `.tbl` file)
- Verify you specified the correct tool with `-t` (hmmscan vs hmmsearch)
- Ensure the FASTA file used matches the one from the HMMER search

### Import errors
Make sure the package is properly installed:
```bash
pip install -e .
```

Or ensure BioPython is available:
```bash
pip install biopython
```
