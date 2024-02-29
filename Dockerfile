# Set the base image to debian jessie
FROM ubuntu:latest
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update
RUN apt-get -y install wget
RUN apt-get -y install perl
RUN apt-get -y install python3
RUN apt-get -y install pip
RUN pip install biopython pandas numpy

#get the binary, open and install? - need to figure out how to install silently. 
RUN apt-get -y install genometools
RUN apt-get -y install bedtools