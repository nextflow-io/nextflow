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

version ?= none

build: dist/docker/amd64
	cp ../nextflow .
	docker buildx build --platform linux/amd64 --output=type=docker --progress=plain --tag nextflow/nextflow:${version} --build-arg TARGETPLATFORM=linux/amd64 .

build-arm: dist/docker/arm64
	cp ../nextflow .
	docker buildx build --platform linux/arm64 --output=type=docker --progress=plain --tag nextflow/nextflow:${version} --build-arg TARGETPLATFORM=linux/arm64 .

release: build
	docker push nextflow/nextflow:${version}
	#
	docker tag nextflow/nextflow:${version} public.cr.seqera.io/nextflow/nextflow:${version}
	docker push public.cr.seqera.io/nextflow/nextflow:${version}

#Static builds can now be found at:
#
#  Linux:   https://download.docker.com/linux/static
#  MacOS:   https://download.docker.com/mac/static
#  Windows: https://download.docker.com/win/static

dist/docker/amd64:
	mkdir -p dist/linux/amd64
	curl -fsSL https://download.docker.com/linux/static/stable/x86_64/docker-17.09.0-ce.tgz | tar --strip-components=1 -xvz -C dist/linux/amd64

dist/docker/arm64:
	mkdir -p dist/linux/arm64
	curl -fsSL https://download.docker.com/linux/static/stable/aarch64/docker-17.09.0-ce.tgz | tar --strip-components=1 -xvz -C dist/linux/arm64

