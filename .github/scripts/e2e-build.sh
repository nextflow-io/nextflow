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
set -e

# Assembles the Nextflow launcher docker build context (Dockerfile + runtime)
# from a checkout at $REPO_ROOT. Deliberately touches no secrets: $REPO_ROOT may
# be an untrusted PR checkout, and this script itself must always be invoked
# from a trusted checkout (see e2e.yml) so a PR can't rewrite what runs here.
: "${REPO_ROOT:?REPO_ROOT must point at the checkout to build}"
CTX="$REPO_ROOT/.github/test-e2e"

rm -rf "$CTX/.nextflow" && mkdir -p "$CTX/.nextflow"
(cd "$REPO_ROOT"
export NXF_PLUGINS_DIR=$PWD/build/plugins
make releaseInfo assemble installScratch
)

cp -r "$REPO_ROOT/build/plugins" "$CTX/.nextflow/"
cp "$REPO_ROOT/nextflow" "$CTX/" && chmod +x "$CTX/nextflow"
cp "$REPO_ROOT/modules/nextflow/src/main/resources/META-INF/build-info.properties" "$CTX/"
source "$CTX/build-info.properties"

if [ -z "$version" ]; then
    echo "Error: version is empty or missing"; exit 1
fi
if [ -z "$build" ]; then
    echo "Error: build is empty or missing"; exit 1
fi
if [ -z "$commitId" ]; then
    echo "Error: commitId is empty or missing"; exit 1
fi

echo "version  : $version"
echo "build    : $build"
echo "commit id: $commitId"
