#!/bin/bash
#
# See 
# https://circleci.com/docs/nightly-builds/
# 
set -e
set -u
set -o pipefail

NXF_VER=${1:-}

function trigger_build() {

	local PROJECT=$1
	local BRANCH=${2:-master}
	echo "Triggering CircleCI build ${PROJECT} [${BRANCH}] -- Nextflow ${NXF_VER:-latest}"

	local trigger_build_url=https://circleci.com/api/v1/project/${PROJECT}/tree/${BRANCH}?circle-token=${CIRCLE_TOKEN}

	local params=''
	[[ $NXF_VER ]] && params+="\"NXF_VER\": \"$NXF_VER\","
	params+="\"RUN_NIGHTLY_BUILD\": \"true\","
	params+="\"FUNCTIONAL_TEST_TARGET\": \"staging-dawn-435.herokuapp.com\""
	
	local post_data=$(cat <<EOF
	{
	  "build_parameters": { $params }
	}
	EOF)

	curl \
	--header "Accept: application/json" \
	--header "Content-Type: application/json" \
	--data "${post_data}" \
	--request POST ${trigger_build_url}

	sleep 1 
}


trigger_build nextflow-io/tests
trigger_build nextflow-io/examples
trigger_build nextflow-io/hello
trigger_build nextflow-io/rnatoy
trigger_build nextflow-io/rnaseq-nf
trigger_build cbcrg/piper-nf
trigger_build cbcrg/kallisto-nf
trigger_build cbcrg/lncRNA-Annotation-nf
trigger_build cbcrg/mta-nf
trigger_build cbcrg/ampa-nf
trigger_build cbcrg/shootstrap
trigger_build CRG-CNAG/CalliNGS-NF
#trigger_build cbcrg/grape-nf -- NO TESTS



