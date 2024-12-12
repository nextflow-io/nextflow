#!/usr/bin/env bash
set -e

# -----------------------------------------------------------------------------
# Nextflow release script
#
# This is the orchestration script for Nextflow releases.
#
# It is intended to be run by 'headless' CI environments (eg Github Actions)
# to execute a release.You probably don't want to run this script directly.
#
# Instead, use the `make release` command to be guided through the process.
# -----------------------------------------------------------------------------

# set defaults for unspecified env vars
export GH_ORG=${GH_ORG:-'nextflow-io'}
export GH_USER=${GH_USER:-'pditommaso'} # TODO - use a service user for releases
export GH_USER_EMAIL=${GH_USER_EMAIL:-'paolo.ditommaso@gmail.com'}

export MAVEN_PUBLISH_URL=${MAVEN_PUBLISH_URL:-'s3://maven.seqera.io'}
export PLUGINS_INDEX_JSON=${PLUGINS_INDEX_JSON:-'https://github.com/nextflow-io/plugins/main/plugins.json'}
export S3_RELEASE_BUCKET=${S3_RELEASE_BUCKET:-'www2.nextflow.io'}
export SEQERA_CONTAINER_REGISTRY=${SEQERA_CONTAINER_REGISTRY:-'public.cr.seqera.io'}

# change to the project root
cd "$(dirname "$0")/.."

# build artifacts
make distribution

# tag release
./release/tag-release.sh

# deploy to maven
./release/deploy-to-maven.sh

# deploy to S3
./release/deploy-to-s3.sh

# deploy to docker
./release/deploy-to-docker.sh

# deploy to github
./release/deploy-to-github.sh

# deploy plugins
./release/deploy-plugins-to-maven.sh
./release/deploy-plugins-to-github.sh
./release/update-plugins-index.sh

# finally, publish the new launcher
./release/publish-launcher-script.sh
