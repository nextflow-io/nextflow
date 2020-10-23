#!/bin/bash
TEST_JDK=${TEST_JDK:=8}
X_BRANCH=${TRAVIS_BRANCH:-${CIRCLE_BRANCH:-'master'}}
X_PULL_REQUEST=${TRAVIS_PULL_REQUEST:-false}
## handle github variables
[[ $GITHUB_REF ]] && X_BRANCH=$(echo $GITHUB_REF | awk '{n=split($1,A,"/"); print A[n]}')
[[ $GITHUB_EVENT_NAME == pull_request ]] && X_PULL_REQUEST=true

if [[ $X_PULL_REQUEST != false ]] && [ ${TEST_JDK} -gt 8 ]; then
  echo "Skipping integration tests on PR and JDK>8"
  exit 0
fi

if [[ $X_BRANCH != master && $X_BRANCH != testing ]] && [ ${TEST_JDK:=8} -gt 8 ]; then
  echo "Integration tests are only executed on branch 'master' or 'testing' for JDK>8"
  exit 0
fi

export WITH_DOCKER='-with-docker'
export NXF_CMD=$PWD/nextflow;
export CAPSULE_LOG=none
export TEST_JDK
unset JAVA_TOOL_OPTIONS # this variable mess-up Capsule loader Java version parsing
(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
