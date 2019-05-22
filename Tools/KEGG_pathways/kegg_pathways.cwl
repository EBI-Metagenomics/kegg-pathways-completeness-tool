#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: kegg:latest

baseCommand: ['python', '/give_pathways.py']

inputs:
  input_table:
    type: File
    inputBinding:
      separate: true
      prefix: -i

stdout: $(inputs.input_table.nameroot)_pathways.txt

outputs:
  output_pathways:
    type: stdout