import numpy as np
import os
import argparse
import sys
import pdb

def parsing(input_file, outdir):
    dict_contigs = {}

    # reading all annotations
    with open(input_file, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            contig = line[3]
            kegg_annotation = line[0]
            if contig not in dict_contigs:
                dict_contigs[contig] = []
            dict_contigs[contig].append(kegg_annotation)

    # leave unique records and save
    path_output = os.path.join(outdir, os.path.basename(input_file).split('.')[0]+'_parsed.txt')
    with open(path_output, 'w+') as file_out:
        for key in dict_contigs:
            #dict_contigs[key] = np.unique(dict_contigs[key])
            file_out.write('\t'.join([key]+list(dict_contigs[key]))+'\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generates file with KEGG orthologs for each contig")
    parser.add_argument("-i", "--input", dest = "input_file", help = "Tab deliminated file with hmmscan results", required = True)
    parser.add_argument("-o", "--outdir", dest = "outdir", help = "Relative path to directory where you want the output file to be stored (default: cwd)", default = ".")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        parsing(args.input_file, args.outdir)