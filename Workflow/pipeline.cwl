#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  input_fasta:
    type: File

outputs:
  hmmscan_out:
    outputSource: hmmscan/output_table
    type: File
  modification_out:
    outputSource: tab_modification/output_with_tabs
    type: File
  parsing_hmmscan_out:
    outputSource: parsing_hmmscan/output_table
    type: File
  kegg_pathways_out:
    outputSource: kegg_pathways/output_pathways
    type: File


steps:
  hmmscan:
    in:
      seqfile: input_fasta
    out:
      - output_table
      - stdout
      - stderr
    run: ../Tools/Hmmscan/hmmscan.cwl

  tab_modification:
    in:
      input_table: hmmscan/output_table
    out:
      - output_with_tabs
    run: ../Tools/Modification/modification_table.cwl

  parsing_hmmscan:
    in:
      table: tab_modification/output_with_tabs
    out:
      - output_table
      - stdout
      - stderr
    run: ../Tools/Parsing_hmmscan/parsing_hmmscan.cwl

  kegg_pathways:
    in:
      input_table: parsing_hmmscan/output_table
    out:
      - output_pathways
  run: ../Tools/KEGG_pathways/kegg_pathways.cwl