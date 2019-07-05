FROM continuumio/miniconda3

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

RUN conda install prokka && conda clean -a
