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

ifdef module 
mm = :${module}:
else 
mm = 
endif 

clean:
	rm -rf .nextflow*
	rm -rf work 
	rm -rf modules/nextflow/.nextflow*
	rm -rf modules/nextflow/work
	rm -rf build
	rm -rf buildSrc/build
	rm -rf modules/*/build
	rm -rf plugins/*/build
	./gradlew clean

compile:
	./gradlew compile exportClasspath
	@echo "DONE `date`"

assemble:
	./gradlew buildInfo compile assemble

check:
	./gradlew check

#
# install compiled artifacts in Maven local dir
# 
install:
	BUILD_PACK=1 ./gradlew installLauncher publishToMavenLocal -Dmaven.repo.local=${HOME}/.nextflow/capsule/deps/

#
# Show dependencies try `make deps config=runtime`, `make deps config=google`
#
deps:
	./gradlew -q ${mm}dependencies --configuration ${config}

deps-all:
	./gradlew -q dependencyInsight --configuration ${config} --dependency ${module}

#
# Refresh SNAPSHOTs dependencies
#
refresh:
	./gradlew --refresh-dependencies 

#
# Run all tests or selected ones
#
test:
ifndef class
	./gradlew ${mm}test
else
	./gradlew ${mm}test --tests ${class}
endif

#
# Run smoke tests
#
smoke:
	NXF_SMOKE=1 ./gradlew ${mm}test

#
# Upload JAR artifacts to Maven Central
#
upload:
	./gradlew upload

#
# Create self-contained distribution package
#
pack:
	BUILD_PACK=1 ./gradlew pack

#
# Upload NF launcher to nextflow.io web site
#
deploy:
	BUILD_PACK=1 ./gradlew deploy

#
# Close artifacts uploaded to Maven central
#
close:
	./gradlew closeAndReleaseRepository

#
# Upload final package to GitHub
#
release:
	BUILD_PACK=1 ./gradlew release

#
# Create and upload docker image distribution
#
dockerImage:
	BUILD_PACK=1 ./gradlew dockerImage

#
# Create local docker image
#
dockerPack:
	BUILD_PACK=1 ./gradlew publishToMavenLocal dockerPack -Dmaven.repo.local=${PWD}/build/docker/.nextflow/capsule/deps/


upload-plugins:
	./gradlew plugins:upload

publish-index:
	./gradlew plugins:publishIndex
