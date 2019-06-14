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

  parsing_hmmscan_out:
    outputSource: parsing_hmmscan/output_table
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
