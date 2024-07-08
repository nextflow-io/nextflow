FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Ben Sherman <bentshermann@gmail.com>"

ARG MAKE_VERSION=4.3
ARG PYTHON_VERSION=3.12
ARG GIT_VERSION=2.45.2

# Set the working directory
WORKDIR /app

# Install required packages and clean up
RUN micromamba install --yes --name base --channel conda-forge \
    make=${MAKE_VERSION} \
    python=${PYTHON_VERSION} \
    git=${GIT_VERSION} && \
    micromamba clean --all --yes

# Copy requirements.txt to the working directory
COPY requirements.txt .

# Install Python packages from requirements.txt
RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate && \
    pip install -r requirements.txt && \
    micromamba clean --all --yes
