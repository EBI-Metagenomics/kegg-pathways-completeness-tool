#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


label: "Biosequence analysis using profile hidden Markov models"

requirements:
  DockerRequirement:
    dockerPull: hmmscan_kegg:latest
  InlineJavascriptRequirement: {}

baseCommand: ["hmmscan"]

arguments:

  #- prefix: -E
  #  valueFrom: "0.001"
  #  position: 2
  - valueFrom: --noali
    position: 1

  - prefix: --domtblout
    valueFrom: $(inputs.seqfile.nameroot)_hmmscan.tbl
    position: 2

  - valueFrom: /db/merged
    position: 3

inputs:

  seqfile:
    type: File
    inputBinding:
      position: 4
      separate: true

stdout: stdout.txt
stderr: stderr.txt

outputs:
  stdout: stdout
  stderr: stderr
  output_table:
    type: File
    outputBinding:
      glob: "*hmmscan.tbl"