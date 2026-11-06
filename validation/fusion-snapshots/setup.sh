#!/bin/bash
#
# Copyright 2013-2026, Seqera Labs
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Creates the FIS role, the experiment template and the Batch job role. Usage: setup.sh <bucket>
# Use a dedicated spot queue and do not run the test concurrently: the experiment
# interrupts every running spot instance tagged fis-target=nf-snapshot-e2e.
set -euo pipefail
cd "$(dirname "$0")"
export BUCKET=${1:?usage: setup.sh <bucket>}
# FIS templates are regional: must be the Batch queue region, not the profile default
export AWS_REGION=${AWS_REGION:?set AWS_REGION to the Batch queue region}
export ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

aws iam create-role --role-name nf-fis-snapshot-e2e \
    --assume-role-policy-document "$(envsubst '$ACCOUNT_ID $AWS_REGION' < fis-role-trust.json)" >/dev/null
aws iam put-role-policy --role-name nf-fis-snapshot-e2e --policy-name fis --policy-document file://fis-role-policy.json
sleep 10  # IAM propagation, FIS rejects a role it cannot see yet
export TEMPLATE_ID=$(aws fis create-experiment-template \
    --cli-input-json "$(envsubst '$ACCOUNT_ID' < fis-template.json)" \
    --query experimentTemplate.id --output text)

JOB_ROLE=$(aws iam create-role --role-name nf-job-snapshot-e2e --assume-role-policy-document file://job-role-trust.json \
    --query Role.Arn --output text)
aws iam put-role-policy --role-name nf-job-snapshot-e2e --policy-name e2e \
    --policy-document "$(envsubst '$TEMPLATE_ID $BUCKET $ACCOUNT_ID $AWS_REGION' < job-role-policy.json)"

echo "export E2E_FIS_TEMPLATE=$TEMPLATE_ID E2E_JOB_ROLE=$JOB_ROLE"
