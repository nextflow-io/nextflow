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
#   train.sh --verify <bundle home> <work dir>    run the training pipeline with the cache of the bundle
#
set -euo pipefail

verify=''
[[ $1 == --verify ]] && verify=1 && shift
home=$(cd "$1" && pwd -P)
if [[ $verify ]]; then
    work=$2
    [[ -f $home/lib/nextflow.aot ]] || { echo "Missing AOT cache: $home/lib/nextflow.aot" >&2; exit 1; }
else
    output=$2
    work=$3
    [[ ! -e $home/lib/nextflow.aot ]] || { echo "The bundle already has an AOT cache: $home/lib/nextflow.aot" >&2; exit 1; }
fi
training=$(cd "$(dirname "$0")/training" && pwd -P)

mkdir -p "$work"
work=$(mktemp -d "$(cd "$work" && pwd -P)/run.XXXXXX")
cp -R "$training/." "$work/"
cd "$work"

# isolate the run from the Nextflow settings of the environment, except the JVM options for --verify
jvm_args=${NXF_JVM_ARGS:-}
for var in $(compgen -e); do [[ $var == NXF_* ]] && unset "$var"; done
export NXF_HOME=$work/home
export NXF_WORK=$work/work
export NXF_CACHE_DIR=$work
export NXF_OFFLINE=true
export NXF_ANSI_LOG=false
export NXF_DISABLE_CHECK_LATEST=true
if [[ $verify ]]; then
    NXF_JVM_ARGS=$jvm_args "$home/bin/nextflow" info
    NXF_JVM_ARGS=$jvm_args "$home/bin/nextflow" run main.nf
else
    mkdir -p "$(dirname "$output")"
    # classes that can't be archived (e.g. generated or in an error state) are reported as warnings
    NXF_JVM_ARGS="-XX:AOTCacheOutput=$output.tmp -Xlog:aot=error" "$home/bin/nextflow" run main.nf
    mv -f "$output.tmp" "$output"
fi
