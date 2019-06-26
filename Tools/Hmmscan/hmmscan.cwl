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

  - valueFrom: "--cut_ga"
    position: 2
  - valueFrom: --noali
    position: 1

  - prefix: --domtblout
    valueFrom: $(inputs.seqfile.nameroot)_hmmscan.tbl
    position: 3


inputs:

  seqfile:
    type: File
    inputBinding:
      position: 5
      separate: true
  data:
    type: Directory?
    default:
      class: Directory
      path:  ../../Files/db/
      location: ../../Files/db/
      listing: []
      basename: db
    inputBinding:
      valueFrom: $(self.path)/db_kofam.hmm  # $(self.listing[0].dirname)/db_kofam.hmm
      position: 4

stdout: stdout.txt
stderr: stderr.txt

outputs:
  stdout: stdout
  stderr: stderr
  output_table:
    type: File
    outputBinding:
      glob: "*hmmscan.tbl"