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
# Agent examples against the nf-agent-pi runner image built from THIS tree.
#
# The plugin jar names its runner image by the plugin VERSION, which is bumped only at release,
# so between a build-context change and its release the image the jar declares is the PREVIOUS
# build - and the new one exists nowhere a Wave or k8s run could pull it from. Publish it to the
# staging registry under a tag derived from the checksum of its build context instead, and run
# the examples against that. The tag is content-addressed, so `push` is a no-op for a build
# context that an earlier run already published.
#
# Requires OPENAI_API_KEY and a Docker login to the staging registry.
set -e

ROOT=$(cd "$(dirname "$0")/.." && pwd)
STAGE=${NF_AGENT_PI_STAGE:-public.cr.stage-seqera.io/nextflow}

tag=$("$ROOT/plugins/nf-agent-pi/build-image.sh" context-tag)
image=$STAGE/nf-agent-pi:$tag
echo "Agent runner image: $image"

# amd64 only: it is what this runner executes the examples on, and building it natively needs no
# QEMU. The release still builds both architectures
"$ROOT/plugins/nf-agent-pi/build-image.sh" push -r "$STAGE" -t "$tag" -P linux/amd64

# validate.sh runs the development launcher, which needs the exported classpath
make -C "$ROOT" compile

echo "Test agent examples on the local executor"
AGENT_VALIDATION_DIR=$PWD/agent-validation \
  "$ROOT/examples/agents/validate.sh" -m local -r -i "$image" \
    01_structured-output 02_two-agents 03_skills
