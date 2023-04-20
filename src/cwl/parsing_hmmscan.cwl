#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/

label: "Taking KOs corresponding contigs from hmmscan result"

hints:
  DockerRequirement:
    dockerPull: kegg_analysis:latest

requirements:
  InlineJavascriptRequirement: {}

baseCommand: ['parsing_hmmscan.py']

inputs:
  table:
    format: edam:format_3475
    type: File
    inputBinding:
      separate: true
      prefix: -i
  fasta:
    type: File
    inputBinding:
      separate: true
      prefix: -f

stdout: stdout.txt
stderr: stderr.txt

outputs:

  output_table:
    type: File
    format: edam:format_3475 # TXT
    outputBinding:
      glob: "*_parsed*"


$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-http.rdf
's:author': 'Ekaterina Sakharova'
's:copyrightHolder': EMBL - European Bioinformatics Institute
's:license': "https://www.apache.org/licenses/LICENSE-2.0"