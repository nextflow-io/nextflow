FROM mambaorg/micromamba:1.3.1

MAINTAINER Ben Sherman <bentshermann@gmail.com>

RUN micromamba install --yes --name base --channel conda-forge \
      make=4.3 \
      python=3.7 \
      conda-forge:git=2.45.0 \
      && \
    micromamba clean --all --yes

COPY requirements.txt .

RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate && \
    pip install -r requirements.txt
