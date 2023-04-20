#!/bin/bash

wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/graphs.pkl.gz \
   ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_class.txt.gz \
   ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_names.txt.gz

gunzip graphs.pkl.gz all_pathways_class.txt.gz all_pathways_names.txt.gz