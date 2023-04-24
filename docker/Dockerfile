FROM alpine:3.14

LABEL Maintainer="Ekaterina Sakharova <kates@ebi.ac.uk>"

WORKDIR pwd
COPY requirements.txt /

# Dependances
RUN apk add --no-cache --update-cache --update python3 gzip && \
apk add --no-cache --update py3-pip && \
apk add --no-cache --update bash && \
pip3 install --upgrade pip && \
pip3 install -r /requirements.txt

# ADD KEGG files to generate pathways
RUN mkdir pathways_data
COPY pathways_data/all_pathways_class.txt graphs/graphs.pkl pathways_data/

# Add .py scripts
RUN mkdir /tools
COPY bin/give_pathways.py docker/run_pathways.sh /tools/

RUN chmod -R a+rwx /tools
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/tools