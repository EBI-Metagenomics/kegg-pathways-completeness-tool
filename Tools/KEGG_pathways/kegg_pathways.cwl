#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: kegg:latest

baseCommand: ['python', '/give_pathways.py']
arguments: ['-o', '/']

inputs:
  input_table:
    type: File
    inputBinding:
      separate: true
      prefix: -i

stdout: stdout.txt

outputs:
  output_pathways_summary:
    type: File
    outputBinding:
      glob: summary_pathways.txt

  output_pathways_matching:
    type: File
    outputBinding:
      glob: matching_ko_pathways.txt

  output_pathways_missing:
    type: File
    outputBinding:
      glob: missing_ko_pathways.txt

  stdout: stdout

  out_folder:
    type: Directory
    outputBinding:
      glob: "Contigs"