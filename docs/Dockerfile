FROM mambaorg/micromamba:1.3.1

MAINTAINER Ben Sherman <bentshermann@gmail.com>

RUN micromamba install --yes --name base --channel conda-forge \
      make=4.3 \
      python=3.7 \
      && \
    micromamba clean --all --yes

RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate && \
    pip install sphinx==5.3.0 myst-parser==0.18.1 sphinx-rtd-theme==1.1.1
