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
  union_by_contigs:
    outputSource: union_by_contigs/output_table
    type: File
  kegg_pathways_summary:
    outputSource: kegg_pathways/output_pathways_summary
    type: File
  kegg_pathways_matching:
    outputSource: kegg_pathways/output_pathways_matching
    type: File
  kegg_pathways_missing:
    outputSource: kegg_pathways/output_pathways_missing
    type: File
  kegg_contigs:
    outputSource: kegg_pathways/out_folder
    type: Directory

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

  union_by_contigs:
    in:
      table: parsing_hmmscan/output_table
    out:
      - output_table
      - stdout
      - stderr
    run: ../Tools/Union_by_contigs/union_by_contigs.cwl

  kegg_pathways:
    in:
      input_table: union_by_contigs/output_table
    out:
      - output_pathways_summary
      - output_pathways_matching
      - output_pathways_missing
      - out_folder
      - stdout
    run: ../Tools/KEGG_pathways/kegg_pathways.cwl