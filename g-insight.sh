#!/bin/bash
CONF=${2:-'compile'}
set -x
./gradlew -q dependencyInsight --configuration compile --dependency $1
