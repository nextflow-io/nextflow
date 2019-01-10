#
#  Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

version ?= none

build: dist/docker 
	cp ../nextflow .
	docker build -t nextflow/nextflow:${version} .

dist/docker:
	curl -fsSLO https://get.docker.com/builds/Linux/x86_64/docker-17.03.1-ce.tgz
	mkdir dist
	tar --strip-components=1 -xvzf docker-17.03.1-ce.tgz -C dist
	rm docker-17.03.1-ce.tgz

release: build
	docker tag nextflow/nextflow:${version} nextflow/nextflow:latest
	docker push nextflow/nextflow:${version}
	docker push nextflow/nextflow:latest
