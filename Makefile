#
#  Copyright 2013-2024, Seqera Labs
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.


#
# add parameters to the command line as `param=value` for example:
# make deps config=runtime
# 

config ?= compileClasspath

gradle = ./gradlew --console=plain

ifdef module 
mm = :${module}:
else 
mm = 
endif 

compile:
	$(gradle) compile exportClasspath
	@echo "DONE `date`"

assemble:
	$(gradle) buildInfo compile assemble

releaseInfo:
	$(gradle) releaseInfo

check:
	$(gradle) check

clean:
	rm -rf .nextflow*
	rm -rf work 
	rm -rf modules/nextflow/.nextflow*
	rm -rf modules/nextflow/work
	rm -rf build
	rm -rf modules/*/build
	rm -rf plugins/*/build
	$(gradle) clean

#
# install compiled artifacts in Maven local dir
# 
install:
	BUILD_PACK=1 \
	$(gradle) installLauncher publishToMavenLocal installPlugin

installScratch:
	BUILD_PACK=1 \
	$(gradle) installScratch publishToMavenLocal installPlugin

#
# Show dependencies try `make deps config=runtime`, `make deps config=google`
#
deps:
	$(gradle) -q ${mm}dependencies --configuration ${config}

deps-all:
	$(gradle) -q dependencyInsight --configuration ${config} --dependency ${module}

#
# Refresh SNAPSHOTs dependencies
#
refresh:
	$(gradle) --refresh-dependencies 

#
# Run all tests or selected ones
#
test:
ifndef class
	$(gradle) ${mm}test
else
	$(gradle) ${mm}test --tests ${class}
endif

#
# Run smoke tests
#
smoke:
	NXF_SMOKE=1 $(gradle) ${mm}test

#
# Upload JAR artifacts to Maven Central
#
upload:
	$(gradle) upload

#
# Create self-contained distribution package
#
pack:
	BUILD_PACK=1 $(gradle) pack

#
# Upload NF launcher to nextflow.io web site
#
deploy:
	BUILD_PACK=1 $(gradle) deploy

#
# Close artifacts uploaded to Maven central
#
close:
	$(gradle) closeAndReleaseRepository

#
# Upload final package to GitHub
#
release:
	BUILD_PACK=1 $(gradle) release

#
# Create and upload docker image distribution
#
dockerImage:
	BUILD_PACK=1 $(gradle) dockerImage

#
# Create local docker image
#
dockerPack:
	BUILD_PACK=1 $(gradle) publishToMavenLocal dockerPack -Dmaven.repo.local=${PWD}/build/docker/.nextflow/capsule/deps/ installPlugin

release-plugins:
	$(gradle) releasePluginToRegistryIfNotExists

#
# Publish the `pi` agent runner image, which is the distribution unit of the nf-agent-pi
# runtime. Runs first in release.sh, because it is the only release step that reaches a
# third-party registry and it must not be able to abort a release that has already
# published something. A no-op when the tag is already published.
#
release-agent-image:
	$(gradle) :plugins:nf-agent-pi:releaseImageIfNotExists

publish-artifacts:
	$(gradle) publishAllPublicationsToSeqeraRepository

