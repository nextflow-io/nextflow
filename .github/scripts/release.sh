#!/usr/bin/env bash
set -e

# build artifacts
make distribution

# tag release
./tag-release.sh

# deploy to maven
./deploy-to-maven.sh

# deploy to S3
./deploy-to-s3.sh

# deploy to docker
./deploy-to-docker.sh

# deploy to github
./deploy-to-github.sh

# deploy plugins
./deploy-plugins-to-maven.sh
./deploy-plugins-to-github.sh
./update-plugins-index.sh

# finally, publish the distribution
./publish-release.sh
