#
#  Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
#
#  This file is part of Nextflow.
#
#  Nextflow is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Nextflow is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.


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

dag:
	@./gradlew -q generateDependencyGraph
	@echo "Check out the DAG at this file: dependency-graph.png"