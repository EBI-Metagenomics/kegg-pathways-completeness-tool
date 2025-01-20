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
    Script to convert hmm-table to tab-separated
"""

import sys
import argparse
from ..utils import __version__

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Convert hmm-table to tab-separated")
    parser.add_argument("-i", "--input", dest="input", help="Input file", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Output file", required=True)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args(argv)


class HmmerTableModifier:
    def __init__(
        self,
        input_file: str,
        output_file: str,
    ):
        """
        Output of hmm-table has columns separated with different number of spaces (to keep humanreadable format)
        :param input_file: output of hmm-table
        :param output_file: name of tab-separated table
        """
        self.input_file = input_file
        self.output_file = output_file

    def modify_table(self):
        with open(self.input_file, 'r') as file_in, open(self.output_file, 'w') as file_out:
            for line in file_in:
                if line.startswith('#'):
                    continue
                line = list(filter(None, line.strip().split(' ')))
                modified_line = '\t'.join(line[:22] + [' '.join(line[22:])])
                file_out.write(modified_line + '\n')


def main():
    args = parse_args(sys.argv[1:])
    hmmtable_table_modifier = HmmerTableModifier(
        input_file=args.input,
        output_file=args.output,
    )
    hmmtable_table_modifier.modify_table()


if __name__ == "__main__":
    main()
