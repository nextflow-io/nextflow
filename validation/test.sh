#!/bin/bash 
set -e 
trap "exit" INT

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}

#
# Tests
#
(
  cd ../tests/checks; 
  bash run.sh
)

#
# Hello 
#
git clone https://github.com/nextflow-io/hello
( 
  cd hello; 
  $NXF_CMD run .
  $NXF_CMD run . -resume
)

#
# Rna-Toy
#
git clone https://github.com/nextflow-io/rnatoy
(
  cd rnatoy; 
  $NXF_CMD run . -with-docker 
  $NXF_CMD run . -with-docker -resume 
)

#
# AMPA-NF
#
git clone https://github.com/cbcrg/ampa-nf
docker pull cbcrg/ampa-nf
(
  cd ampa-nf; 
  $NXF_CMD run . -with-docker 
  $NXF_CMD run . -with-docker -resume 
)

#
# MTA-NF
#
git clone https://github.com/cbcrg/mta-nf
docker pull cbcrg/mta-nf
(
  cd mta-nf;
  $NXF_CMD run . -with-docker --seq tutorial/small.fa --ntree 5 --msa clustalw
  $NXF_CMD run . -with-docker --seq tutorial/small.fa --ntree 5 --msa clustalw -resume
)

#
# GRAPE-NF
#
git clone  https://github.com/cbcrg/grape-nf
docker pull cbcrg/grape-nf
(
  cd grape-nf; 
  $NXF_CMD run . -with-docker 
  $NXF_CMD run . -with-docker -resume 
)

#
# PIPER-NF
#
git clone https://github.com/cbcrg/piper-nf
docker pull cbcrg/piper-nf
(
  cd piper-nf; 
  $NXF_CMD run . -with-docker 
  $NXF_CMD run . -with-docker -resume 
)



