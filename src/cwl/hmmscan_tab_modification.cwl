cwlVersion: v1.0
class: CommandLineTool
$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/

requirements:
  ResourceRequirement:
    coresMax: 1
    ramMin: 100  # just a default, could be lowered

hints:
  DockerRequirement:
    dockerPull: mgnify/pipeline-v5.python3:latest

inputs:
  input_table:
    format: edam:format_3475
    type: File
    inputBinding:
      prefix: '-i'

baseCommand: [ hmmscan_tab.py ]  # old was with sed

arguments:
  - valueFrom: $(inputs.input_table.nameroot).tsv
    prefix: -o

outputs:
  output_with_tabs:
    type: File
    format: edam:format_3475  # TXT
    outputBinding:
      glob: "*.tsv"


$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf
  - http://edamontology.org/EDAM_1.16.owl

's:author': 'Ekaterina Sakharova'
's:copyrightHolder': EMBL - European Bioinformatics Institute
's:license': "https://www.apache.org/licenses/LICENSE-2.0"