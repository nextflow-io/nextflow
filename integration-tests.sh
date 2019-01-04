#!/bin/bash

if [[ $TRAVIS_PULL_REQUEST == true ]] && [ ${TEST_JDK:=8} -gt 8 ]; then
  echo "Skipping integration tests on PR and JDK>8"
  exit 0
fi

if [[ $TEST_JDK ]]; then
    curl -o install-jdk.sh https://raw.githubusercontent.com/sormuras/bach/master/install-jdk.sh
    chmod +x install-jdk.sh
    source ./install-jdk.sh --cacerts --feature $TEST_JDK
fi

export WITH_DOCKER='-with-docker'
export NXF_CMD=$PWD/nextflow; 
(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
