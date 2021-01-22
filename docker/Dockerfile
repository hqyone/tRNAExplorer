FROM ubuntu:18.04

#####################################################################################
############################# How to use the docker image ###########################
#####################################################################################
## Using docker file to build image
# cd <docker folder>
## Build the image
# docker build -t <trna> .  
## Run the image <trna> with an interactive way
# docker run --name <trna_cont> -it <trna>  
## Run tRNAExplorer or tRNA_db_maker in the container
# root@ff8c666f8aca:/tools# tRNAExplorer.py
# root@ff8c666f8aca:/tools# tRNA_db_maker.py  

## Run container for jupyter notebook (visulization)
# docker run -p 7777:8888 -v $(pwd):/ext -it trna
## lauch the notebook for visualization analysis processed data
# root@ff8c666f8aca:/tools# jupyter noteboook
## connect server using the URL : <IP>:7777, example:http://192.168.0.10:7777
## use "trna" as password to login.


######################## Docker image for tRNAExplorer ##############################
#####################################################################################

LABEL MAINTAINER Quanyuan He <hqyone@gmail.com>

#####################################################################################
###### Environment settings and dependency installation - DO NO EDIT ################
#####################################################################################

USER root
ENV TERM=vt100
ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies including python, wget, x11, R, etc
RUN apt-get update && \
    apt-get -y --no-install-recommends install apt-utils build-essential dialog gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libcurl4-gnutls-dev && \   
    apt-get -y --no-install-recommends install libssl-dev libxml2-dev xorg-dev wget less ca-certificates openjdk-8-jre-headless python3.7 unzip git vim

RUN apt-get -y --no-install-recommends install python3-pip python3-setuptools samtools trnascan-se ncbi-blast+ trimmomatic
#RUN ln -s /usr/bin/python3.7 /usr/bin/python3
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade setuptools && \
    pip3 install --no-cache-dir setuptools jupyter scipy numpy pandas matplotlib pyfaidx && \
    pip3 install --no-cache-dir seaborn pysam pathlib scipy

RUN pip3 install plotly==4.14.3

#####################################################################################
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
RUN git clone https://github.com/hqyone/tRNAExplorer.git /tools/tRNAExplorer && \
    printf 't_path="/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"\nblastn="/usr/bin/blastn"\nmkdb="/usr/bin/makeblastdb"'>/tools/tRNAExplorer/init && \
    printf 'out_dir=/tools/tRNAExplorer/test/output\nsample_tsv=/tools/tRNAExplorer/test/samples\ntrna_anno_bed=/tools/tRNAExplorer/test/trna_db/hg38_tRNA_60.bed"'>/tools/tRNAExplorer/test/output/visual_config.tsv && \
    chmod +x /tools/tRNAExplorer/tRNAExplorer.py && \
    chmod +x /tools/tRNAExplorer/tRNA_db_maker.py

## setting jupyter notebook passwords and ports 
RUN mkdir /root/.jupyter/
## use "trna" as password
RUN printf "c.NotebookApp.allow_root = True\nc.NotebookApp.ip = '*'\nc.NotebookApp.notebook_dir = '/tools/tRNAExplorer/visualization'\nc.NotebookApp.password = u'argon2:\$argon2id\$v=19\$m=10240,t=10,p=8\$aC5yGldsL/djCLPvF6KlaQ\$RQGgWufIGWcTdHY5FbKtAw'">/root/.jupyter/jupyter_notebook_config.py
ENV PATH="/tools/tRNAExplorer:${PATH}"


## Downloads human and mouse genome files

## WORKDIR /data/hg38
# RUN wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P /data/hg38/
# RUN unzip /data/hg38/hg38.fa.gz -d /src/data/hg38
# RUN samtools faidx /data/hg38/hg38.fa

## WORKDIR /src/data/mm10
# RUN wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -P /data/mm10/
# RUN unzip /data/mm10/mm10.fa.gz -d /data/mm10
# RUN samtools faidx /data/mm10/mm10.fa


