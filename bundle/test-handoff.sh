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
# Tests the hand-off of the installed launcher (the `nextflow` script of this repository) to a bundle:
# serves the bundle from a file:// mirror with the layout of the release server, lets the launcher
# install it into a fresh NXF_HOME, and runs the training pipeline through it with the AOT cache
# (-XX:AOTMode=on: the JVM fails to start if the cache is rejected).
#
# Usage:
#   test-handoff.sh <bundle tar.gz> <work dir> [<command>...]
#
# The optional command is run after the training pipeline, from the current directory, in the
# environment of the hand-off and with NXF_CMD set to the launcher, e.g. to run integration tests.
#
set -euo pipefail

function die() { echo "$*" >&2; exit 1; }

cwd=$PWD
tarball=$(cd "$(dirname "$1")" && pwd -P)/$(basename "$1")
root=$(cd "$(dirname "$0")/.." && pwd -P)
# nextflow-<version>-<os>-<arch>
name=$(basename "$tarball" .tar.gz)
version=${name#nextflow-}
version=${version%-*-*}
platform=${name#nextflow-"$version"-}
[[ -f $tarball.sha256 ]] || die "Missing checksum file: $tarball.sha256"

mkdir -p "$2"
work=$(mktemp -d "$(cd "$2" && pwd -P)/handoff.XXXXXX")
mkdir -p "$work/mirror/v$version" "$work/run"
cp "$tarball" "$tarball.sha256" "$work/mirror/v$version/"
cp -R "$root/bundle/training/." "$work/run/"
cd "$work/run"

# isolate the run from the Nextflow settings of the environment
for var in $(compgen -e); do [[ $var == NXF_* ]] && unset "$var"; done
export NXF_HOME=$work/home
export NXF_WORK=$work/run/work
export NXF_ANSI_LOG=false
export NXF_DISABLE_CHECK_LATEST=true
export NXF_VER=$version
export NXF_JAVA_BUNDLE=true
export NXF_JAVA_BUNDLE_BASE=file://$work/mirror
export NXF_JVM_ARGS=-XX:AOTMode=on

info=$("$root/nextflow" info -d)
head -n 5 <<< "$info"
installed=$NXF_HOME/framework/$version/nextflow-$version-$platform
[[ -L $installed && -x $installed/runtime/bin/java ]] || die "The launcher did not install the bundle in $installed"
# the AOT cache is only passed to the JVM by the launcher of the bundle, with its runtime
grep -qF -- "-XX:AOTCache=$(cd "$installed" && pwd -P)/lib/nextflow.aot" <<< "$info" \
    || die "The launcher did not hand off to the bundle in $installed"
"$root/nextflow" run main.nf

if [[ $# -gt 2 ]]; then
    cd "$cwd"
    NXF_CMD=$root/nextflow "${@:3}"
fi
