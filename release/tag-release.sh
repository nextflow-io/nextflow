#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/.."

# read the nextflow version
read -r NF_VERSION<VERSION

echo "
----------------------
-- Creating git tag --
----------------------
"

# create & push a git tag
git tag "v$NF_VERSION"
git push origin "v$NF_VERSION"

echo "Done"
