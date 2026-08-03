FROM debian:stable-slim

# Tools required by the test suite: blastp (blast.nf, blast-parallel.nf),
# t_coffee (blast-parallel.nf, watch.nf), AMPA.pl (ampa.nf).
# Math::CDF is not packaged by Debian, so it is the only source build.
# gnuplot is required: AMPA.pl exits 2 if absent, even with -noplot.
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -qq -y --no-install-recommends \
      ncbi-blast+ \
      t-coffee \
      libmath-round-perl \
      gnuplot-nox \
      build-essential \
      cpanminus && \
    cpanm --notest Math::CDF && \
    apt-get purge -qq -y build-essential cpanminus && \
    apt-get autoremove -qq -y && \
    rm -rf /var/lib/apt/lists/* /root/.cpanm

COPY bin/AMPA.pl /usr/local/bin/

CMD ["/bin/bash"]
