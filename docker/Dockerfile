FROM alpine:3.5

LABEL Maintainer="Ekaterina Sakharova <kates@ebi.ac.uk>"

# Dependances
RUN apk add --no-cache --update-cache --update python3 gzip && \
apk add --no-cache --update py3-pip && \
apk add --no-cache --update bash && \
pip3 install --upgrade pip && \
pip3 install networkx

# Downloading KEGG help files to generate pathways
RUN mkdir help_files
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/graphs.pkl.gz -P help_files && gzip -d help_files/graphs.pkl.gz
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_names.txt.gz -P help_files && gzip -d help_files/all_pathways_names.txt.gz
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/all_pathways_class.txt.gz -P help_files && gzip -d help_files/all_pathways_class.txt.gz

# Add .py scripts
RUN mkdir /tools
ADD Tools/give_pathways.py /tools/give_pathways.py
ADD run_pathways.sh /tools/run_pathways.sh

RUN chmod -R a+rwx /tools
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/tools