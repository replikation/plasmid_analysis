FROM replikation/debian_basic
MAINTAINER Christian Brandt <christian.jena@gmail.com>
# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda
RUN rm -f Miniconda3-latest-Linux-x86_64.sh
ENV PATH /conda/bin:$PATH
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install abricate
RUN abricate --check
RUN abricate-get_db --db ncbi --force
RUN abricate-get_db --db plasmidfinder --force
RUN mkdir -p /conda/db/transposon
COPY sequences /conda/db/transposon
RUN abricate --setupdb
