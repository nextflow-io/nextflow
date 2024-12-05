#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/../.."

echo "Publishing nextflow jars to maven"

# the release process should have already built the jars, so to avoid re-compiling everything
# we can tell gradle to skip all non publish/publication related tasks
./gradlew publishToMaven \
  $( ./gradlew publishToMaven --dry-run | grep -iv 'publish\|publication' | awk '/^:/ { print "-x" $1 }')

echo "Done"
