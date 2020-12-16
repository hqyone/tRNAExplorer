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
    apt-get -y --no-install-recommends install apt-utils build-essential dialog gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev libcurl4-gnutls-dev libssl-dev libxml2-dev xorg-dev wget less openjdk-8-jre-headless python3.7 unzip

RUN apt-get -y  --no-install-recommends install r-base
RUN apt-get -y install texinfo texlive-latex-base texlive-fonts-recommended texlive-fonts-extra
#RUN ln -s /usr/bin/python3.7 /usr/bin/python3

#####################################################################################
#################### Install software and download database #########################
#####################################################################################

RUN mkdir /src/
WORKDIR /src/

## Install BLASTN
# RNN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
# RUN tar -xf ncbi-blast-2.11.0+-x64-linux.tar.gz
# RUN ln -s /src/ncbi-blast-2.11.0+/bin/blastn /usr/bin/blastn

RUN sudo apt -y install ncbi-blast+

## Install BLASTN
RNN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
RUN unzip Trimmomatic-0.39.zip
RUN ln -s /src/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/bin/trimmomatic-0.39.jar