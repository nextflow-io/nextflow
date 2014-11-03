#!/usr/bin/env bash
set -e
set -x 

#
# verify Amazon auth codes are availables
#
[[ -z 'AWS_ACCESS_KEY' ]] && { echo 'Missing $AWS_ACCESS_KEY variable'; exit 1; }
[[ -z 'AWS_SECRET_KEY' ]] && { echo 'Missing $AWS_SECRET_KEY variable'; exit 2; }

mkdir -p $HOME/bin
cat <<EOF >> $HOME/.bash_profile
export AWS_ACCESS_KEY=$AWS_ACCESS_KEY
export AWS_SECRET_KEY=$AWS_SECRET_KEY
EOF

#
# Lauch docker and pull the container when DOCKER variable is defined
#
[[ "$DOCKER_IMAGE" ]] && docker pull $DOCKER_IMAGE

# the bucket name
AWS_S3BUCKET=${AWS_S3BUCKET:-'nxf-cluster'}

#
# Install NEXTFLOW and launch it
#
if [[ "$NXF_VER" ]]; then
  version="v$NXF_VER"
else
  version='latest'
fi

curl -fsSL http://www.nextflow.io/releases/${version}/nextflow  > $HOME/nextflow
chmod +x $HOME/nextflow
bash -x $HOME/nextflow node -bg \
  -cluster.join "s3:$AWS_S3BUCKET" \
  -cluster.interface eth0

# save the environment for debugging 
env | sort > .boot.env


