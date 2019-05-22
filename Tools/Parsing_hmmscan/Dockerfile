FROM python:3

LABEL maintainer="Ekaterina Sakharova <kates@ebi.ac.uk>"

RUN pip install pandas
RUN pip install biopython

##############################################################
# Dockerfile Version:   18.06.0-ce
# Software:             parsing_hmmscan.py
# Software Version:     1.0.0
##############################################################
COPY parsing_hmmscan.py /parsing_hmmscan.py


CMD ["python", "/parsing_hmmscan.py"]