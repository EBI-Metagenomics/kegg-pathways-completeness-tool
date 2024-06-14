#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Script to convert hmmscan table to tab-separated
"""

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert fastq to fasta")
    parser.add_argument("-i", "--input", dest="input", help="Input file", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Output file", required=True)

    args = parser.parse_args()

    with open(args.input, 'r') as file_in, open(args.output, 'w') as file_out:
        for line in file_in:
            if line.startswith('#'):
                continue
            line = list(filter(None, line.strip().split(' ')))
            modified_line = '\t'.join(line[:22] + [' '.join(line[22:])])
            file_out.write(modified_line + '\n')

if __name__ == "__main__":
    main()
