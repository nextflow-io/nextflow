#!/usr/bin/env bash
set -e
set -x 

#
# verify Amazon auth codes are availables
#
[[ -z 'AWS_ACCESS_KEY' ]] && { echo 'Missing $AWS_ACCESS_KEY variable'; exit 1; }
[[ -z 'AWS_SECRET_KEY' ]] && { echo 'Missing AWS_SECRET_KEY variable'; exit 2; }

#
# Lauch docker and pull the container when DOCKER variable is defined
#
if [ -z "DOCKER" ]; then
    docker pull $DOCKER
fi

# the bucket name
AWS_S3BUCKET=${AWS_S3BUCKET:-'nxf-cluster'}

#
# Install NEXTFLOW and launch it
#
curl -fsSL http://get.nextflow.io > nextflow && chmod +x nextflow
./nextflow -bg \
  -daemon.name gridgain \
  -daemon.interface eth0 \
  -daemon.join s3:$AWS_S3BUCKET

# save the environment for debugging 
env | sort > .boot.env


