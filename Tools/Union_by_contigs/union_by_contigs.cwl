#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


label: "Unite KOs related to the same contig"

requirements:
  DockerRequirement:
    dockerPull: kegg_union_by_contigs:latest
  InlineJavascriptRequirement: {}

baseCommand: ['python', '/union_by_contigs.py']

inputs:
  table:
    type: File
    inputBinding:
      separate: true
      prefix: -i

stdout: stdout.txt
stderr: stderr.txt


outputs:
  stdout: stdout
  stderr: stderr

  output_table:
    type: File
    outputBinding:
      glob: "union_ko_contigs.txt"