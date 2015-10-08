#!/bin/bash
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

MODULE=${2:-''}
CONFIG=${1:-'compile'}

if [[ $1 == 'help' || $1 == '-h' ]]; then
  echo Show dependencies for a specific configuration and sub-project
  echo USAGE: g-deps [configuration] [sub-project]
  echo 
  echo Examples:
  echo ' deps.sh runtime                  --  Show runtime dependencies in main project'
  echo ' deps.sh testCompile nxf-dnanexus --  Show test dependencies in nxf-dnanexus sub-project' 
  exit 1 
fi

set -x
./gradlew -q $MODULE:dependencies --configuration $CONFIG