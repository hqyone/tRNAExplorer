FROM ubuntu:18.04

######################## Docker image for tRNAExplorer ##############################
#####################################################################################

LABEL MAINTAINER Quanyuan He <hqyone@gmail.com>

#####################################################################################
###### Environment settings and dependency installation - DO NO EDIT ################
#####################################################################################

ENV TERM=vt100
ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies including python, wget, x11, R, etc
RUN apt-get update && \
    apt-get -y --no-install-recommends install apt-utils build-essential dialog gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libcurl4-gnutls-dev && \   
    apt-get -y --no-install-recommends install libssl-dev libxml2-dev xorg-dev wget less ca-certificates openjdk-8-jre-headless python3.7 unzip git

RUN apt-get -y --no-install-recommends install python3-pip python3-setuptools samtools trnascan-se ncbi-blast+ trimmomatic
#RUN ln -s /usr/bin/python3.7 /usr/bin/python3
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade setuptools && \
    pip3 install --no-cache-dir setuptools jupyter scipy numpy pandas matplotlib pyfaidx && \
    pip3 install --no-cache-dir seaborn pysam pathlib scipy

RUN pip3 install plotly==4.14.3

#####################################################################################ƒ
#################### Install software and download database #########################
#####################################################################################

RUN mkdir /data
RUN mkdir /data/hg38
RUN mkdir /data/mm10
RUN mkdir /ext
RUN mkdir /tools
WORKDIR /tools

## Install trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -P /tools/
RUN unzip /tools/Trimmomatic-0.39.zip -d /tools

## Install tRNAExplorer
RUN git clone https://github.com/hqyone/tRNAExplorer.git /tools/tRNAExplorer
RUN printf 't_path="/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"\nblastn="/usr/bin/blastn"\nmkdb="/usr/bin/makeblastdb"'>/tools/tRNAExplorer/init
ENV PATH="/tools/tRNAExplorer:${PATH}"

## Downloads human and mouse genome files

# WORKDIR /src/data/hg38
# RUN wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P /src/data/hg38/
# RUN unzip /src/data/hg38/hg38.fa.gz -d /src/data/hg38
# RUN samtools faidx /src/data/hg38/hg38.fa


# WORKDIR /src/data/mm10
# RUN wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -P /src/data/mm10/
# RUN unzip /src/data/mm10/mm10.fa.gz -d /src/data/mm10
# RUN samtools faidx /src/data/mm10/mm10.fa

## Usage the file
#
# use 'docker build -t trna .' to build the image
# use 'docker run --name trna_cont -it trna python3 /tools/tRNAExplorer' to run the container in an interactive way 
# 


