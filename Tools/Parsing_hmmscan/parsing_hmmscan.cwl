#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


label: "Biosequence analysis using profile hidden Markov models"

requirements:
  DockerRequirement:
    dockerPull: parsing_hmmscan:latest
  InlineJavascriptRequirement: {}

baseCommand: ['python', '/parsing_hmmscan.py']

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
      glob: "*_parsed.txt"