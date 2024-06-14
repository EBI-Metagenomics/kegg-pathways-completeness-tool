#!/usr/bin/env python3
import os
import argparse
import sys
from Bio import SeqIO

def get_dir_contigs(input_fasta):
    dict_contigs = {}
    seq_records = SeqIO.parse(input_fasta, "fasta")
    for line in seq_records:
        if line.name not in dict_contigs:
            dict_contigs[line.name] = []
    print(len(dict_contigs))
    return dict_contigs


def parsing(dict_contigs, input_file, outdir):
    # reading all annotations
    with open(input_file, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            contig = line[0]  #line[3] - hmmscan
            kegg_annotation = line[3]  # line[0] - hmmscan
            contig_in_fasta = [name for name in dict_contigs if name in contig]
            if len(contig_in_fasta) == 0:
                print(contig)
                continue
            elif len(contig_in_fasta) == 1:
                dict_contigs[contig_in_fasta[0]].append(kegg_annotation)
            else:
                print('strange contig')

    # leave unique records and save
    path_output = os.path.join(outdir, os.path.basename(input_file)+'_parsed')
    with open(path_output, 'w+') as file_out:
        for key in dict_contigs:
            if len(dict_contigs[key]) != 0:
                file_out.write('\t'.join([key]+list(dict_contigs[key]))+'\n')

def main():
    parser = argparse.ArgumentParser(description="Generates file with KEGG orthologs for each contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Tab deliminated file with hmmscan results",
                        required=True)
    parser.add_argument("-f", "--fasta", dest="fasta_file", help="Filtered fasta file with initial names of contigs",
                        required=True)
    parser.add_argument("-o", "--outdir", dest="outdir", help="Relative path to directory where you want the output file to be stored (default: cwd)", default = ".")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        parsing(get_dir_contigs(args.fasta_file), args.input_file, args.outdir)


if __name__ == "__main__":
    main()
