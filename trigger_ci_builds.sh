#!/bin/bash
#
# See 
# https://docs.travis-ci.com/user/triggering-builds/
# 
set -e
set -u
set -o pipefail

NXF_VER=${1:-}

function trigger_build() {

	local PROJECT=$1
	local BRANCH=${2:-master}
	echo "Triggering Travis build ${PROJECT} [${BRANCH}] -- Nextflow ${NXF_VER:-latest}"

	local body='{ "request": { "branch":"master" '
	if [[ $NXF_VER ]]; then
	  body+=', "config":{ "env": {"NXF_VER":"'$NXF_VER'"}} '
	fi 
	body+='}}'

	curl -s -X POST \
	-H "Content-Type: application/json" \
	-H "Accept: application/json" \
	-H "Travis-API-Version: 3" \
	-H "Authorization: token $TRAVIS_TOKEN" \
	-d "$body" \
	https://api.travis-ci.org/repo/$(echo -n $PROJECT | sed 's@/@%2F@g')/requests

	sleep 1 
}


trigger_build nextflow-io/tests
trigger_build nextflow-io/examples
trigger_build nextflow-io/hello
trigger_build nextflow-io/rnatoy
trigger_build nextflow-io/rnaseq-nf
trigger_build cbcrg/piper-nf
trigger_build cbcrg/kallisto-nf
trigger_build cbcrg/mta-nf
trigger_build cbcrg/ampa-nf
trigger_build cbcrg/shootstrap
trigger_build CRG-CNAG/CalliNGS-NF
