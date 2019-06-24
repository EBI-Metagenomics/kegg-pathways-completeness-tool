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

  - prefix: --domE
    valueFrom: "0.00001"
    position: 2
  - valueFrom: --noali
    position: 1

  - prefix: --domtblout
    valueFrom: $(inputs.seqfile.nameroot)_hmmscan.tbl
    position: 3

  - valueFrom: /db/merged
    position: 4

inputs:

  seqfile:
    type: File
    inputBinding:
      position: 5
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