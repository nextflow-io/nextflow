#
#  Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

config ?= compile

ifdef module 
mm = :${module}:
else 
mm = 
endif 

compile:
	./gradlew compile exportClasspath
	@echo "DONE `date`"

clean:
	./gradlew clean

assemble:
	./gradlew compile assemble

check:
	./gradlew check

install:
	./gradlew installLauncher install -Dmaven.repo.local=${HOME}/.nextflow/capsule/deps/ -x signArchives

deps:
	./gradlew -q ${mm}dependencies --configuration ${config}

deps-all:
	./gradlew -q dependencyInsight --configuration ${config} --dependency ${module}

refresh:
	./gradlew --refresh-dependencies 

test:
ifndef class
	./gradlew ${mm}test
else
	./gradlew ${mm}test --tests ${class}
endif

pack:
	./gradlew packAll

deploy:
	./gradlew deploy

close:
	./gradlew closeRepository promoteRepository
	
release:
	./gradlew release	
	
dockerPack:
	./gradlew install dockerPack -Dmaven.repo.local=${PWD}/build/docker/.nextflow/capsule/deps/ -x signArchives
