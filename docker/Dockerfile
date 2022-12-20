FROM --platform=linux/amd64 continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y time
RUN apt-get install -y file
RUN conda install -y -c conda-forge -c bioconda metaplatanus=1.3.1
