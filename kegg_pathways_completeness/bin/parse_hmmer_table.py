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

"""
Script to process HMMER output tables and extract KO annotations per contig.
"""

import os
import tempfile
import click
from Bio import SeqIO
from importlib.metadata import version, PackageNotFoundError


def get_version():
    """Get package version from installed metadata"""
    try:
        return version("kegg-pathways-completeness")
    except PackageNotFoundError:
        return "unknown"


class HmmerTableConverter:
    """Converts HMMER table output to tab-separated format."""

    def __init__(self, input_file: str, output_file: str):
        """
        Output of hmm-table has columns separated with different number of spaces
        (to keep human-readable format). This class converts it to tab-separated format.

        :param input_file: output of hmm-table
        :param output_file: name of tab-separated table
        """
        self.input_file = input_file
        self.output_file = output_file

    def convert(self):
        """Convert space-separated HMMER table to tab-separated format."""
        with open(self.input_file, "r") as file_in, open(
            self.output_file, "w"
        ) as file_out:
            for line in file_in:
                if line.startswith("#"):
                    continue
                line = list(filter(None, line.strip().split(" ")))
                modified_line = "\t".join(line[:22] + [" ".join(line[22:])])
                file_out.write(modified_line + "\n")


class HmmerTableParser:
    """Parses tab-separated HMMER table to generate KOs per contig."""

    def __init__(
        self,
        input_table: str,
        input_fasta: str,
        output_dir: str,
        hmmtool: str,
        output_basename: str = None,
    ):
        """
        Script parses hmmer tab-separated file. It generates a list of KOs per each contig.
        Contig names in hmmer output can be different from initial fasta file.
        In order to keep contig names consistent, it is required to parse initial fasta file too.

        :param input_table: hmmer tab-separated table
        :param input_fasta: fasta file that was used in hmmer command
        :param output_dir: name of output directory
        :param hmmtool: hmmscan / hmmsearch
        :param output_basename: optional custom basename for output file (if not provided, uses input_table basename)
        """
        self.input_table = input_table
        self.input_fasta = input_fasta
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.hmmtool = hmmtool
        self.output_basename = output_basename

    def get_dir_contigs(self):
        """Extract contig names from FASTA file."""
        dict_contigs = {}
        seq_records = SeqIO.parse(self.input_fasta, "fasta")
        for line in seq_records:
            if line.name not in dict_contigs:
                dict_contigs[line.name] = []
        click.echo(f"Found {len(dict_contigs)} contigs in FASTA file")
        return dict_contigs

    def choose_columns(self, line):
        """Select appropriate columns based on HMMER tool used."""
        # return contig, KO
        if self.hmmtool == "hmmsearch":
            return line[0], line[3]
        elif self.hmmtool == "hmmscan":
            return line[3], line[0]
        else:
            raise ValueError(f"Incorrect HMM tool specified: {self.hmmtool}")

    def parse(self, dict_contigs):
        """Parse the tab-separated HMMER table and assign KOs to contigs."""
        # reading all annotations
        with open(self.input_table, "r") as file_in:
            for line in file_in:
                line = line.strip().split("\t")
                contig, kegg_annotation = self.choose_columns(line)
                contig_in_fasta = [name for name in dict_contigs if name in contig]
                if len(contig_in_fasta) == 0:
                    click.echo(f"Warning: Contig {contig} not found in FASTA file")
                    continue
                elif len(contig_in_fasta) == 1:
                    dict_contigs[contig_in_fasta[0]].append(kegg_annotation)
                else:
                    click.echo(f"Warning: Ambiguous contig match for {contig}")

        # leave unique records and save
        if self.output_basename:
            basename = self.output_basename
        else:
            basename = os.path.basename(self.input_table)
        path_output = os.path.join(self.output_dir, basename + "_parsed")
        with open(path_output, "w+") as file_out:
            for key in dict_contigs:
                if len(dict_contigs[key]) != 0:
                    file_out.write("\t".join([key] + list(dict_contigs[key])) + "\n")
        click.echo(f"Output saved to: {path_output}")
        return path_output


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True),
    help="Input HMMER domtblout file",
)
@click.option(
    "-f",
    "--fasta",
    "fasta_file",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file with protein sequences used in HMMER search",
)
@click.option(
    "-t",
    "--tool",
    "hmm_tool",
    required=True,
    type=click.Choice(["hmmscan", "hmmsearch"], case_sensitive=False),
    help="HMMER tool used (hmmscan or hmmsearch)",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    required=True,
    type=click.Path(),
    help="Output file for KO annotations per contig",
)
@click.option(
    "--save-intermediate",
    "intermediate_file",
    type=click.Path(),
    help="Save intermediate tab-separated table to this file (optional)",
)
@click.version_option(version=get_version(), prog_name="parse_hmmer_table")
def main(input_file, fasta_file, hmm_tool, output_file, intermediate_file):
    """
    Process HMMER domtblout output to extract KO annotations per contig.

    This tool converts HMMER's domain table output format to a tab-separated
    format and then parses it to generate a file with KO annotations for each
    contig. The output format is: contig_name<TAB>KO1<TAB>KO2<TAB>...

    \b
    Example:
    parse_hmmer_table -i hmmer_output.tbl -f proteins.faa -t hmmscan -o result.tsv

    \b
    With intermediate file:
    parse_hmmer_table -i hmmer.tbl -f proteins.faa -t hmmscan \\
        -o result.tsv --save-intermediate hmmer_intermediate.tsv
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file) or "."
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Convert to tab-separated format
    click.echo("Step 1: Converting HMMER table to tab-separated format...")

    if intermediate_file:
        # User wants to save the intermediate file
        temp_tsv = intermediate_file
        keep_temp = True
    else:
        # Create temporary file in same directory as output
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, dir=output_dir
        ) as tmp_file:
            temp_tsv = tmp_file.name
        keep_temp = False

    try:
        converter = HmmerTableConverter(input_file=input_file, output_file=temp_tsv)
        converter.convert()
        if intermediate_file:
            click.echo(f"Tab-separated table saved to: {temp_tsv}")
        else:
            click.echo("Tab-separated table created (temporary)")

        # Step 2: Parse the tab-separated table
        click.echo("\nStep 2: Parsing tab-separated table to extract KOs per contig...")
        contigs = {}
        seq_records = SeqIO.parse(fasta_file, "fasta")
        for line in seq_records:
            if line.name not in contigs:
                contigs[line.name] = []
        click.echo(f"Found {len(contigs)} contigs in FASTA file")

        # Reading all annotations
        with open(temp_tsv, "r") as file_in:
            for line in file_in:
                line = line.strip().split("\t")
                # Choose columns based on tool
                if hmm_tool == "hmmsearch":
                    contig, kegg_annotation = line[0], line[3]
                elif hmm_tool == "hmmscan":
                    contig, kegg_annotation = line[3], line[0]
                else:
                    raise ValueError(f"Incorrect HMM tool specified: {hmm_tool}")

                contig_in_fasta = [name for name in contigs if name in contig]
                if len(contig_in_fasta) == 0:
                    click.echo(f"Warning: Contig {contig} not found in FASTA file")
                    continue
                elif len(contig_in_fasta) == 1:
                    contigs[contig_in_fasta[0]].append(kegg_annotation)
                else:
                    click.echo(f"Warning: Ambiguous contig match for {contig}")

        # Save output
        with open(output_file, "w+") as file_out:
            for key in contigs:
                if len(contigs[key]) != 0:
                    file_out.write("\t".join([key] + list(contigs[key])) + "\n")

        click.echo("\nProcessing complete!")
        click.echo(f"Final output: {output_file}")
    finally:
        # Clean up temporary file only if user didn't request to save it
        if not keep_temp and os.path.exists(temp_tsv):
            os.remove(temp_tsv)


if __name__ == "__main__":
    main()
