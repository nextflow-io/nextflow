FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Ben Sherman <bentshermann@gmail.com>"

ARG MAKE_VERSION=4.3
# Note netlify latest currently supported python version is 3.8
ARG PYTHON_VERSION=3.8
ARG GIT_VERSION=2.45.2

RUN micromamba install --yes --name base --channel conda-forge \
    make=${MAKE_VERSION} \
    python=${PYTHON_VERSION} \
    conda-forge:git=${GIT_VERSION} && \
    micromamba clean --all --yes

COPY requirements.txt .

RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate && \
    pip install -r requirements.txt && \
    micromamba clean --all --yes
