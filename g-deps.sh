#!/bin/bash
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