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

# cleanup
rm -rf .nextflow && mkdir .nextflow
# copy nextflow dependencies
(cd ..
export NXF_PLUGINS_DIR=$PWD/build/plugins
make releaseInfo assemble installScratch
)

# copy nextflow plugins
cp -r ../build/plugins .nextflow/
# copy nextflow launcher script
cp ../nextflow . && chmod +x nextflow
cp ../modules/nextflow/src/main/resources/META-INF/build-info.properties .
source build-info.properties

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
base=${base:-'public.cr.seqera.io/platform/nf-launcher:j21-base'}
base_dind=${base_dind:-'public.cr.seqera.io/platform/nf-launcher:j21-dind-base'}
repository=${repository:-'public.cr.seqera.io/snapshots/nextflow-scratch'}
image=${repository}:${tag}

docker buildx build \
  --platform linux/amd64 \
  --push \
  --progress=plain \
  --tag ${image} \
  --build-arg TARGETPLATFORM=linux/amd64 \
  .
echo "Nextflow snapshots launcher image $image"

#
# Create ephemeral containers with the scratch image and the base Platform launcher images.
# As in Platform, the dind launcher (including the Docker client) is used by cloud compute
# environments, the plain launcher by all others
#
launcher=$(wave -i ${base} --include ${image} --platform linux/amd64 --config-env NXF_HOME=/.nextflow)
launcher_dind=$(wave -i ${base_dind} --include ${image} --platform linux/amd64 --config-env NXF_HOME=/.nextflow)
echo "Running Platform tests using image launcher: $launcher"
echo "Running Platform tests using image launcher (dind): $launcher_dind"

# determine the e2e test environment: use $E2E_ENVIRONMENT if set (workflow_dispatch),
# otherwise check $COMMIT_MESSAGE set by the GitHub action workflow, falling back
# to the commit message in the git repo
if [ -n "$E2E_ENVIRONMENT" ]; then
  ENVIRONMENT="$E2E_ENVIRONMENT"
  echo "Environment [from input]: $ENVIRONMENT"
else
  if [ -z "$COMMIT_MESSAGE" ]; then
    COMMIT_MESSAGE=$(git show -s --format='%s')
    echo "Commit message [from git]: $COMMIT_MESSAGE"
  else
    echo "Commit message [from gha]: $COMMIT_MESSAGE"
  fi
  if echo "$COMMIT_MESSAGE" | grep -q "\[e2e prod\]"; then
    ENVIRONMENT="production"
  else
    ENVIRONMENT="staging"
  fi
fi

#
# Finally launch the showcase automation
# see https://github.com/seqeralabs/showcase-automation/
showcase_repo=seqeralabs/showcase-automation
workflow_file=seqera-showcase-${ENVIRONMENT}.yml
echo "Launching ${workflow_file}"
dispatch_time=$(date -u +%Y-%m-%dT%H:%M:%SZ)
gh workflow run \
  "${workflow_file}" \
  --repo "${showcase_repo}" \
  -f launch_container=${launcher} \
  -f launch_container_dind=${launcher_dind}

# `gh workflow run` does not return the dispatched run, so look for the newest
# run created since the dispatch, falling back to the workflow runs page
run_url=""
for _ in $(seq 1 12); do
  sleep 5
  run_url=$(gh run list --repo "${showcase_repo}" --workflow "${workflow_file}" --event workflow_dispatch --limit 10 \
    --json url,createdAt --jq "[.[] | select(.createdAt >= \"${dispatch_time}\")] | sort_by(.createdAt) | last | .url // empty" 2>/dev/null || true)
  [ -n "$run_url" ] && break
done
[ -z "$run_url" ] && run_url="https://github.com/${showcase_repo}/actions/workflows/${workflow_file}"
echo "Showcase run: $run_url"

if [ -n "$GITHUB_STEP_SUMMARY" ]; then
  {
    echo "### e2e showcase launched against \`${ENVIRONMENT}\`"
    echo "- Run: ${run_url}"
    echo "- Launcher: \`${launcher}\`"
    echo "- Launcher (dind): \`${launcher_dind}\`"
    echo ""
    echo "CI only launches the test. Open the run to check whether it passed."
  } >> "$GITHUB_STEP_SUMMARY"
fi

