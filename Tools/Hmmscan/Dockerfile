FROM centos:centos7

LABEL maintainer="Ekaterina Sakharova <kates@ebi.ac.uk>"

RUN yum -y update; yum clean all
RUN yum -y install gcc git make wget; yum clean all
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && tar -zxvf hmmer.tar.gz && cd hmmer-* && ./configure && make && make install
#COPY db/ /db

CMD ["hmmscan"]
