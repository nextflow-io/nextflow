#!/usr/bin/env bash
set -e

# change to the project root
cd "$(dirname "$0")/../.."

# read the nextflow version
read -r NF_VERSION<VERSION

# build the docker image
cd docker
make build version="$NF_VERSION"

# push to docker hub
IMAGE_TAG="nextflow/nextflow:$NF_VERSION"
docker push "$IMAGE_TAG"

# push to seqera registry
if [ -z "$SEQERA_REGISTRY" ]; then
  echo "Missing env var SEQERA_REGISTRY"
  exit 1
fi
SEQERA_IMAGE_TAG="$SEQERA_REGISTRY/$IMAGE_TAG"

docker tag "$IMAGE_TAG" "$SEQERA_IMAGE_TAG"
docker push "$SEQERA_IMAGE_TAG"

echo "Done"
