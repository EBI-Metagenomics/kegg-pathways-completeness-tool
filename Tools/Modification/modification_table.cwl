#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: sed_docker:latest

baseCommand: ['sed', '/^#/d; s/ \+/\t/g']

inputs:
  input_table:
    type: File
    inputBinding:
      separate: true
      position: 2

stdout: $(inputs.input_table.nameroot)_tab.tbl

outputs:
  output_with_tabs:
    type: stdout