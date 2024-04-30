FROM mambaorg/micromamba:1.5.8

RUN micromamba install --yes --name base --channel conda-forge \
      make=4.3 \
      python=3.7 \
      git=2.44.0 \
      && \
    micromamba clean --all --yes

COPY requirements.txt .

RUN eval "$(micromamba shell hook --shell=bash)" && \
    micromamba activate && \
    pip install -r requirements.txt
