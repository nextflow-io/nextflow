#!/bin/bash
#
# setup env for integration tests
#
TEST_JDK=${TEST_JDK:=11}
X_BRANCH=${TRAVIS_BRANCH:-${CIRCLE_BRANCH:-'master'}}
X_PULL_REQUEST=${TRAVIS_PULL_REQUEST:-false}
## handle github variables
[[ $GITHUB_REF ]] && X_BRANCH=$(echo $GITHUB_REF | awk '{n=split($1,A,"/"); print A[n]}')
[[ $GITHUB_EVENT_NAME == pull_request ]] && X_PULL_REQUEST=true

if [ "$TEST_JDK" -ge 19 ]; then
  export NXF_ENABLE_VIRTUAL_THREADS=true
fi

export WITH_DOCKER='-with-docker'
export NXF_PLUGINS_DIR=$PWD/build/plugins
export NXF_CMD=$PWD/nextflow;
export CAPSULE_LOG=none
export TEST_JDK
export TEST_MODE

unset JAVA_TOOL_OPTIONS # this variable mess-up Capsule loader Java version parsing
(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
