#!/bin/bash 
set -e 

get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

export NXF_IGNORE_WARN_DSL2=true
export NXF_CMD=${NXF_CMD:-$(get_abs_filename ../launch.sh)}
export TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST:=false}

#
# Tests
#
(
  cd ../tests/
  sudo bash cleanup.sh
  cd checks
  bash run.sh
)

if [[ $TEST_SMOKE == true ]]; then
  echo Skipping tests since TEST_SMOKE flag is true
  exit 0
fi

# disable ansi log to make log more readable
export NXF_ANSI_LOG=false

#
# Hello 
#
git clone https://github.com/nextflow-io/hello
( 
  cd hello; 
  $NXF_CMD run .
  $NXF_CMD run . -resume
)

if [[ $TRAVIS_PULL_REQUEST != false ]]; then
echo Skipping tests requiring secret vars
exit 0
fi

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
# RNASEQ-NF
#
echo nextflow-io/rnaseq-nf
[[ $TOWER_ACCESS_TOKEN ]] && OPTS='-with-tower' || OPTS=''
$NXF_CMD run nextflow-io/rnaseq-nf -with-docker $OPTS
$NXF_CMD run nextflow-io/rnaseq-nf -with-docker $OPTS -resume

#
# AWS Batch tests
#
echo aws batch tests
bash awsbatch.sh

#
# Google Life Sciences
#
if [[ $GOOGLE_SECRET ]]; then
  bash gls.sh
else
  echo "Google Life Science test skipped because GOOGLE_SECRET env var is missing"
fi
