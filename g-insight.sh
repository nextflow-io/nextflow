#!/bin/bash
CONF=${2:-'compile'}
set -x
gradle -q dependencyInsight --configuration compile --dependency $1
