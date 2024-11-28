#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/../.."

# read the nextflow version
read -r NF_VERSION<VERSION

echo "Publishing nextflow release to github"

# create a github (pre)release and attach launcher and dist files
# use --verify-tag to fail if tag doesn't exist
gh release create \
  --prerelease \
  --title "Version $NF_VERSION" \
  --verify-tag \
  "v$NF_VERSION" \
  nextflow \
  "build/releases/nextflow-$NF_VERSION-dist"

echo "Done"
