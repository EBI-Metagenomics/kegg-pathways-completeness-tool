#!/bin/bash

while getopts i: option; do
	case "${option}" in
		i) INPUT=${OPTARG};;
	esac
done

export OUTPUT_PARSING=parsing_table
export OUTPUT_UNION=union_ko_contigs.txt
export OUTDIR=pathways
mkdir -p ${OUTDIR}
export OUTPUT=${OUTDIR}/result

echo "pathways" && \
python3 /tools/give_pathways.py -i ${INPUT} -g /pathways_data/graphs.pkl -c /pathways_data/all_pathways_class.txt -n /pathways_data/all_pathways_names.txt -o ${OUTPUT} && \
echo "finish"