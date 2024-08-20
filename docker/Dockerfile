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

RUN pip install https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool/archive/refs/tags/1.0.5.zip

ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin: