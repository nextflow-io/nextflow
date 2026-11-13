#!/usr/bin/env bash
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
# Creates the AOT cache of a bundle by running the training pipeline with the bundle launcher,
# so that the cache is created with the same JVM options as when it is used.
#
# Usage:
#   train.sh <bundle home> <output cache> <work dir>
#
set -euo pipefail

home=$(cd "$1" && pwd -P)
output=$2
work=$3
[[ ! -e $home/lib/nextflow.aot ]] || { echo "The bundle already has an AOT cache: $home/lib/nextflow.aot" >&2; exit 1; }
training=$(cd "$(dirname "$0")/training" && pwd -P)

mkdir -p "$work"
work=$(mktemp -d "$(cd "$work" && pwd -P)/run.XXXXXX")
cp -R "$training/." "$work/"
cd "$work"

# isolate the run from the Nextflow settings of the environment
for var in $(compgen -e); do [[ $var == NXF_* ]] && unset "$var"; done
export NXF_HOME=$work/home
export NXF_WORK=$work/work
export NXF_CACHE_DIR=$work
export NXF_OFFLINE=true
export NXF_ANSI_LOG=false
export NXF_DISABLE_CHECK_LATEST=true
mkdir -p "$(dirname "$output")"
# classes that can't be archived (e.g. generated or in an error state) are reported as warnings
NXF_JVM_ARGS="-XX:AOTCacheOutput=$output.tmp -Xlog:aot=error" "$home/bin/nextflow" run main.nf -profile test
mv -f "$output.tmp" "$output"
