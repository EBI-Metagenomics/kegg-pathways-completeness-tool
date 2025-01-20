#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import argparse
import sys
from Bio import SeqIO
from ..utils import __version__


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generates file with KEGG orthologs for each contig")
    parser.add_argument("-i", "--input", dest="input_file", help="Tab deliminated file with hmm-table results",
                        required=True)
    parser.add_argument("-f", "--fasta", dest="fasta_file", help="Filtered fasta file with initial names of contigs",
                        required=True)
    parser.add_argument("-o", "--outdir", dest="outdir", default=".",
                        help="Relative path to directory where you want the output file to be stored (default: cwd)")
    parser.add_argument("-t", "--tool", dest="hmm_tool", help="hmmer main command name",
                        choices=['hmmscan', 'hmmsearch'], required=True)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args(argv)


class HmmerTableParser:
    def __init__(
        self,
        input_table: str,
        input_fasta: str,
        output_dir: str,
        hmmtool: str,
    ):
        """
        Script parses hmmer tab-separated file. It generates a list of KOs per each contig.
        Contig names in hmmer output can be different from initial fasta file.
        In order to keep contig names consistent, it is required to parse initial fasta file too.
        :param input_table: hmmer tab-separated table
        :param input_fasta: fasta file that was used in hmmer command
        :param output_dir: name of output directory
        :param hmmtool: hmmscan / hmmsearch
        """
        self.input_table = input_table
        self.input_fasta = input_fasta
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.hmmtool = hmmtool

    def get_dir_contigs(self):
        dict_contigs = {}
        seq_records = SeqIO.parse(self.input_fasta, "fasta")
        for line in seq_records:
            if line.name not in dict_contigs:
                dict_contigs[line.name] = []
        print(len(dict_contigs))
        return dict_contigs

    def choose_columns(self, line):
        # return contig, KO
        if self.hmmtool == 'hmmsearch':
            return line[0], line[3]
        elif self.hmmtool == 'hmmscan':
            return line[3], line[0]
        else:
            raise ValueError(f"Incorrect HMM tool specified: {self.hmmtool}")

    def parsing(self, dict_contigs):
        # reading all annotations
        with open(self.input_table, 'r') as file_in:
            for line in file_in:
                line = line.strip().split('\t')
                contig, kegg_annotation = self.choose_columns(line)
                contig_in_fasta = [name for name in dict_contigs if name in contig]
                if len(contig_in_fasta) == 0:
                    print(contig)
                    continue
                elif len(contig_in_fasta) == 1:
                    dict_contigs[contig_in_fasta[0]].append(kegg_annotation)
                else:
                    print('strange contig')

        # leave unique records and save
        path_output = os.path.join(self.output_dir, os.path.basename(self.input_table)+'_parsed')
        with open(path_output, 'w+') as file_out:
            for key in dict_contigs:
                if len(dict_contigs[key]) != 0:
                    file_out.write('\t'.join([key]+list(dict_contigs[key]))+'\n')


def main():
    args = parse_args(sys.argv[1:])
    hmmer_table_parser = HmmerTableParser(
        input_table=args.input_file,
        input_fasta=args.fasta_file,
        output_dir=args.outdir,
        hmmtool=args.hmm_tool
    )
    contigs = hmmer_table_parser.get_dir_contigs()
    hmmer_table_parser.parsing(contigs)


if __name__ == "__main__":
    main()
