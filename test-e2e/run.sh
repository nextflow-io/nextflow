#
# Copyright 2013-2024, Seqera Labs
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

# cleanup
rm -rf .nextflow && mkdir .nextflow
# copy nextflow dependencies
(cd ..
./gradlew compile assemble
BUILD_PACK=1 ./gradlew installScratch publishToMavenLocal
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

#
# build a scratch container image with assembled newxtflow runtime and plugins
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
  .
echo "Nextflow snapshots launcher image $image"

#
# Create an ephemeral container with the scratch image and base Platform launcher image
#
launcher=$(wave -i ${base} --include ${image} --config-env NXF_HOME=/.nextflow)
echo "Running Platform tests using image launcher: $launcher"

# determining the e2e test environment checking the $COMMIT_MESSAGE
# that is set by GitHub action workflow. If it does not exist fallback
# to the commit message in the git rpeo
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

#
# Finally launch the showcase automation
# see https://github.com/seqeralabs/showcase-automation/
echo "Launching seqera-showcase-${ENVIRONMENT}"
gh workflow run \
  seqera-showcase-${ENVIRONMENT}.yml \
  --repo seqeralabs/showcase-automation \
  -f launch_container=${launcher}

