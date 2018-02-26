#
#  Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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
