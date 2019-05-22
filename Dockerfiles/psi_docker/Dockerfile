FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive
MAINTAINER Christian Brandt <christian.jena@gmail.com>

RUN apt-get update && apt-get install -y \
        git \
        g++ \
        make \
        ncbi-blast+ \
        zlib1g-dev

RUN git clone https://github.com/weizhongli/cdhit.git && \
    cd cdhit && \
    make

ENV PATH /cdhit/psi-cd-hit:$PATH
ENV PATH /cdhit:$PATH
