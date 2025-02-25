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

export WITH_DOCKER='-with-docker'
export NXF_PLUGINS_DIR=$PWD/build/plugins
export NXF_CMD=$PWD/nextflow;
export TEST_JDK
export TEST_MODE

(
 $NXF_CMD info
 cd validation
 bash -x test.sh
)
