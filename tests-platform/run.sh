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

TAG=${version}-${commitId}
base=${base:-'public.cr.seqera.io/platform/nf-launcher:j17-base'}
repository=${repository:-'public.cr.seqera.io/snapshots/nextflow-scratch'}
image=${repository}:${TAG}

docker buildx build \
  --platform linux/amd64 \
  --push \
  --progress=plain \
  --tag ${image} \
  --build-arg TARGETPLATFORM=linux/amd64 \
  .

launcher=$(wave -i ${base} --include ${image} --config-env NXF_HOME=/.nextflow)

echo "Running Platform tests using image launcher: $launcher"

gh workflow run \
  seqera-showcase-staging.yml \
  --repo seqeralabs/showcase-automation \
  -f launch_container=$launcher

#tw launch hello \
#  -w seqeralabs/showcase \
#  --launch-container $launcher
