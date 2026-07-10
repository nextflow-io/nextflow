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

# Build a snapshot Nextflow launcher container and trigger the Seqera Platform
# showcase e2e test, surfacing the dispatched run URL via step outputs/summary.
# ENVIRONMENT ('production'|'staging') is set by the workflow; if unset it falls
# back to grepping the commit message for '[e2e prod]'.
SHOWCASE_REPO=${SHOWCASE_REPO:-'seqeralabs/showcase-automation'}
REPO_ROOT="$(git rev-parse --show-toplevel)"
# docker build context: Dockerfile + assembled runtime (installScratch targets its .nextflow/)
CTX="$REPO_ROOT/.github/test-e2e"

# cleanup
rm -rf "$CTX/.nextflow" && mkdir -p "$CTX/.nextflow"
# copy nextflow dependencies
(cd "$REPO_ROOT"
export NXF_PLUGINS_DIR=$PWD/build/plugins
make releaseInfo assemble installScratch
)

# copy nextflow plugins
cp -r "$REPO_ROOT/build/plugins" "$CTX/.nextflow/"
# copy nextflow launcher script
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

#
# build a scratch container image with assembled nextflow runtime and plugins
#
tag=${version}-${commitId}
base=${base:-'public.cr.seqera.io/platform/nf-launcher:j17-base'}
repository=${repository:-'public.cr.seqera.io/snapshots/nextflow-scratch'}
image=${repository}:${tag}

docker buildx build \
  --platform linux/amd64 \
  --push \
  --progress=plain \
  --tag ${image} \
  --build-arg TARGETPLATFORM=linux/amd64 \
  "$CTX"
echo "Nextflow snapshots launcher image $image"

#
# Create an ephemeral container with the scratch image and base Platform launcher image
#
launcher=$(wave -i ${base} --include ${image} --platform linux/amd64 --config-env NXF_HOME=/.nextflow --config-env NXF_SYNTAX_PARSER=v1)
echo "Running Platform tests using image launcher: $launcher"

# determine the e2e environment: prefer $ENVIRONMENT, else fall back to the commit message
if [ -z "$ENVIRONMENT" ]; then
  [ -z "$COMMIT_MESSAGE" ] && COMMIT_MESSAGE=$(git show -s --format='%s')
  if echo "$COMMIT_MESSAGE" | grep -q "\[e2e prod\]"; then
    ENVIRONMENT="production"
  else
    ENVIRONMENT="staging"
  fi
fi

#
# Launch the showcase automation
# see https://github.com/seqeralabs/showcase-automation/
#
workflow_file="seqera-showcase-${ENVIRONMENT}.yml"
echo "Launching ${workflow_file} in ${SHOWCASE_REPO}"
dispatch_time=$(date -u +%Y-%m-%dT%H:%M:%SZ)
gh workflow run "${workflow_file}" --repo "${SHOWCASE_REPO}" -f launch_container=${launcher}

# `gh workflow run` does not return the dispatched run, so poll for the newest
# run created at/after the dispatch time (falling back to the workflow runs page)
run_url=""
for attempt in $(seq 1 12); do
  sleep 5
  run_url=$(gh run list --repo "${SHOWCASE_REPO}" --workflow "${workflow_file}" --event workflow_dispatch --limit 10 \
    --json url,createdAt --jq "[.[] | select(.createdAt >= \"${dispatch_time}\")] | sort_by(.createdAt) | last | .url // empty" 2>/dev/null || true)
  [ -n "$run_url" ] && break
done
[ -z "$run_url" ] && run_url="https://github.com/${SHOWCASE_REPO}/actions/workflows/${workflow_file}"
echo "Showcase run: $run_url"

# expose results to the calling workflow (step outputs) and the run page (summary)
if [ -n "$GITHUB_OUTPUT" ]; then
  {
    echo "environment=${ENVIRONMENT}"
    echo "launcher_image=${launcher}"
    echo "showcase_url=${run_url}"
  } >> "$GITHUB_OUTPUT"
fi
if [ -n "$GITHUB_STEP_SUMMARY" ]; then
  {
    echo "### e2e showcase launched against \`${ENVIRONMENT}\`"
    echo "- Run: [${run_url}](${run_url})"
    echo "- Launcher image: \`${launcher}\`"
    echo "> CI only launches the test — open the run above to confirm it passed."
  } >> "$GITHUB_STEP_SUMMARY"
fi
