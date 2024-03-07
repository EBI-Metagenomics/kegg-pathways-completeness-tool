#!/bin/bash

export INPUT=false
export LIST_KO=false

while getopts i:l: option; do
    case "${option}" in
        i) INPUT=${OPTARG};;
        l) LIST_KO=${OPTARG};;
    esac
done

export OUTPUT_PARSING=parsing_table
export OUTPUT_UNION=union_ko_contigs.txt
export OUTDIR=pathways
mkdir -p ${OUTDIR}
export OUTPUT=${OUTDIR}/result

if [ $INPUT ]; then
    echo "Processing input file"
    python3 /tools/give_pathways.py -i ${INPUT} -g /pathways_data/graphs.pkl -c /pathways_data/all_pathways_class.txt -n /pathways_data/all_pathways_names.txt -o ${OUTPUT_PARSING}
    echo "Finished processing"
fi

if [ $LIST_KO ]; then
    echo "Processing KO list"
    python3 /tools/give_pathways.py -l ${LIST_KO} -g /pathways_data/graphs.pkl -c /pathways_data/all_pathways_class.txt -n /pathways_data/all_pathways_names.txt -o ${OUTPUT_UNION}
    echo "Finished processing"
fi
