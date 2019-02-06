#!/bin/bash
X_JDK=${TEST_JDK:=8}
X_BRANCH=${TRAVIS_BRANCH:-${CIRCLE_BRANCH:-'master'}}
X_PULL_REQUEST=${TRAVIS_PULL_REQUEST:-false}

if [[ $X_PULL_REQUEST != false ]] && [ ${X_JDK} -gt 8 ]; then
  echo "Skipping integration tests on PR and JDK>8"
  exit 0
fi

if [[ $X_BRANCH != master && $X_BRANCH != testing ]] && [ ${X_JDK:=8} -gt 8 ]; then
  echo "Integration tests are only executed on branch 'master' or 'testing' for JDK>8"
  exit 0
fi

if [ ${X_JDK} -gt 8 ]; then
    curl -o install-jdk.sh https://raw.githubusercontent.com/sormuras/bach/master/install-jdk.sh
    chmod +x install-jdk.sh
    source ./install-jdk.sh --cacerts --feature $X_JDK
fi

export WITH_DOCKER='-with-docker'
export NXF_ANSI_LOG=false
export NXF_CMD=$PWD/nextflow;
(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
