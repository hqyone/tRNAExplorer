FROM ubuntu:18.04

######################## Docker image for tRNAExplorer ##############################
#####################################################################################

MAINTAINER Quanyuan He <hqyone@gmail.com>

#####################################################################################
###### Environment settings and dependency installation - DO NO EDIT ################
#####################################################################################

ENV TERM=vt100
ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies including python, wget, x11, R, etc
RUN apt-get update && \
    apt-get -y --no-install-recommends install apt-utils build-essential dialog gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libcurl4-gnutls-dev libssl-dev libxml2-dev xorg-dev wget less ca-certificates openjdk-8-jre-headless python3.7 unzip
#RUN apt-get -y --no-install-recommends install python-pip samtools trnascan-se ncbi-blast+ trimmomatic
#RUN ln -s /usr/bin/python3.7 /usr/bin/python3
#RUN pip install --update pip && \
#    pip install --no-cache-dir jupyter scipy numpy pandas matplotlib pyfaidx && \
#    pip install --no-cache-dir seaborn pysam HTSeq pathlib plotyly scipy &&

#####################################################################################
#################### Install software and download database #########################
#####################################################################################

RUN mkdir /src/
RUN mkdir /src/data/
RUN mkdir /src/data/hg38
RUN mkdir /src/data/mm10
RUN mkdir /src/tools/
WORKDIR /src/

## Install tRNAExplorer
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip

## Downloads human and mouse genome files
RUN wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P /src/data/hg38/
RUN unzip /src/data/hg38/hg38.fa.gz
RUN samtools faidx /src/data/hg38/hg38.fa
RUN makeblastdb -in /src/data/hg38/hg38.fa -dbtype nucl -title hg38

#S  FSDFS



# WORKDIR /src/data/mm10
# RUN wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -P /src/data/mm10/
# RUN unzip /src/data/mm10/mm10.fa.gz
# RUN samtools faidx /src/data/mm10/mm10.fa
# RUN makeblastdb -in /src/data/mm10/mm10.fa -dbtype nucl -title mm10


## Usage the file
#
# use 'docker build -t ucla .' to build the image
# use 'docker run --name ucla_cont -it ucla /bin/bash' to run the container in an interactive way 

