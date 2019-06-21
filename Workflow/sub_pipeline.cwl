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
  input_table:
    type: File

outputs:

  union_by_contigs:
    outputSource: union_by_contigs/output_table
    type: File


steps:
  tab_modification:
    in:
      input_table: input_table
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

  union_by_contigs:
    in:
      table: parsing_hmmscan/output_table
    out:
      - output_table
      - stdout
      - stderr
    run: ../Tools/Union_by_contigs/union_by_contigs.cwl