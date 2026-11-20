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

#
# Checks that released plugins still link against this Nextflow build: every class, method
# and field they use exists, and no inherited abstract method is left unimplemented.
#
#   plugin-linkage.sh <plugin-id>[@<version>]...
#
# Without a version it checks the latest release. Run from the Nextflow root with Java 25,
# after `./gradlew exportClasspath`.
#
set -euo pipefail

SCRIPTS=$(dirname "$0")
WORK=build/plugin-linkage
REGISTRY=https://registry.nextflow.io/api/v1/plugins
mkdir -p $WORK

[ -s .launch.classpath ] || { echo "Missing .launch.classpath: run ./gradlew exportClasspath"; exit 1; }
core=$(<.launch.classpath)

status=0
for spec in "$@"; do
    id=${spec%@*}
    version=${spec#"$id"}; version=${version#@}
    release=$(curl -sSfL "$REGISTRY/$id" | jq -ce --arg v "$version" '.plugin.releases | if $v == ""
        then map(select(.status == "PUBLISHED")) | max_by(.version | split(".") | map(tonumber? // 0))
        else map(select(.version == $v))[0] end') \
        || { echo "Cannot get release $spec from the registry"; exit 1; }
    version=$(jq -r .version <<<"$release")
    dir=$WORK/$id-$version
    echo "== $id@$version"

    # Registry downloads are public, but need an anonymous token.
    token=$(curl -sSfL "https://auth.cr.seqera.io/auth/token?service=public.cr.seqera.io&scope=repository:nextflow/plugin/$id:pull" | jq -r .token)
    curl -sSfL -H "Authorization: Bearer $token" -o $dir.zip "$(jq -r .url <<<"$release")"
    sha512sum -c --quiet <<<"$(jq -r .sha512sum <<<"$release")  $dir.zip"
    rm -rf $dir && unzip -q $dir.zip -d $dir

    # Classpath: Nextflow core, the parent plugins built from this checkout, and the plugin's own jars.
    # ponytail: parents of parents are not added; no plugin here depends on another one yet.
    cp=$core
    for parent in $(jq -r '.dependsOn[]? | sub("@.*"; "")' <<<"$release"); do
        [ -f plugins/$parent/VERSION ] || { echo "Parent plugin $parent is not in this repo"; exit 1; }
        cp+=":plugins/$parent/build/libs/$parent-$(<plugins/$parent/VERSION).jar"
        cp+=$(printf ':%s' plugins/$parent/build/target/libs/*.jar)
    done
    cp+=$(printf ':%s' $dir/lib/*.jar)

    java $SCRIPTS/Linkage.java "$cp" $dir/classes || status=1
done
exit $status
