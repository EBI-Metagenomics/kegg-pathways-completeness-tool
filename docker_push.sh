#!/usr/bin/env bash
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin

#docker push mgnify/hmmscan_kegg:latest
#docker push mgnify/kegg:latest
#docker push mgnify/parsing_hmmscan:latest
#docker push mgnify/sed_docker:latest
#docker push mgnify/kegg_union_by_contigs:latest