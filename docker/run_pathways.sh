#!/bin/bash

export INPUT=false
export LIST_KO=false
export OUTPUT="pathways_result"

while getopts o:i:l: option; do
    case "${option}" in
        o) OUTPUT=${OPTARG};;
        i) INPUT=${OPTARG};;
        l) LIST_KO=${OPTARG};;
    esac
done

if [ $INPUT ]; then
    echo "Processing input file"
    python3 /tools/give_pathways.py -i ${INPUT} -g /pathways_data/graphs.pkl -c /pathways_data/all_pathways_class.txt -n /pathways_data/all_pathways_names.txt -o ${OUTPUT}
    echo "Finished processing"
fi

if [ $LIST_KO ]; then
    echo "Processing KO list"
    python3 /tools/give_pathways.py -l ${LIST_KO} -g /pathways_data/graphs.pkl -c /pathways_data/all_pathways_class.txt -n /pathways_data/all_pathways_names.txt -o ${OUTPUT}
    echo "Finished processing"
fi
