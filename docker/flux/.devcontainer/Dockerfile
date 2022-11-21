FROM fluxrm/flux-sched:focal

# From project root!
# docker build -f docker/flux/.devcontainer/Dockerfile -t nextflow-flux .
# docker run -it nextflow-flux
# make

LABEL maintainer="Vanessasaurus <@vsoch>"

# Pip not provided in this version
USER root
RUN apt-get update && apt-get install -y python3-venv default-jdk
ENV PYTHONPATH=/usr/lib/flux/python3.8

# For optional usage of Python bindings.
RUN python3 -m pip install IPython
ENV PATH=/env/bin:${PATH}
WORKDIR /code
COPY . /code
ENTRYPOINT /bin/bash
