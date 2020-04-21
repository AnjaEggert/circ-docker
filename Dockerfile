# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

ARG BASE_CONTAINER=jupyter/r-notebook

FROM $BASE_CONTAINER

LABEL maintainer="Anja Eggert <eggert@fbn-dummerstorf.de>"
LABEL description="Circadian analysis"

USER root

RUN apt-get update && \
    apt-get upgrade -y --no-install-recommends \
    apt-get install -y --no-install-recommends make && \
    rm -rf /var/lib/apt/lists/*

# Add channels to install specific R-packages
#RUN conda update -n base conda
RUN conda config --add channels defaults \
 && conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda install --yes r-rocr r-furrr r-rfit r-matrixstats bioconductor-rain bioconductor-qvalue \
 && conda clean --all -f -y \
 && fix-permissions $CONDA_DIR

# Install package from local .tar.gz
RUN Rscript -e 'install.packages("HarmonicRegression_1.9999.tar.gz")'

USER $NB_UID

# Make directory
RUN mkdir /home/jovyan/analysis
