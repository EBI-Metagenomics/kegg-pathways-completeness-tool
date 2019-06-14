FROM python:3
FROM centos:centos7

LABEL maintainer="Ekaterina Sakharova <kates@ebi.ac.uk>"


##### Dependances
RUN yum -y update; yum clean all
RUN yum -y install curl wget fetch perl bzip2 make; yum clean all


##### GNU parallel
RUN (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - pi.dk/3) | bash
RUN echo 'will cite' | parallel --citation &> /dev/null || true


##### Python packages
RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"; python get-pip.py
RUN pip install networkx numpy


##### KEGG pathways downloading
RUN mkdir pathways

# downloading list of names of pathways
RUN curl -s http://rest.kegg.jp/list/module | cut -f1 | cut -c 4- > pathways/list_pathways.txt

# taking DEFINITION file with pathway from KEGG file
# saving in format: <name:pathway>
RUN parallel -k echo -n '{}:' ';' curl -s http://rest.kegg.jp/get/{} '|' grep ^DEFINITION '|' cut -c 13- :::: pathways/list_pathways.txt > pathways/all_pathways.txt
# taking NAME
RUN parallel -k echo -n '{}:' ';' curl -s http://rest.kegg.jp/get/{} '|' grep ^NAME '|' cut -c 13- :::: pathways/list_pathways.txt > pathways/all_pathways_names.txt
# taking CLASS
RUN parallel -k echo -n '{}:' ';' curl -s http://rest.kegg.jp/get/{} '|' grep ^CLASS '|' cut -c 13- :::: pathways/list_pathways.txt > pathways/all_pathways_class.txt


##### Create a Graph for each pathway
ADD make_graphs.py .
RUN python make_graphs.py -i pathways/all_pathways.txt

COPY . /


CMD ["python", "/give_pathways.py"]