#!/usr/bin/env bash
set -e
set -x 

#
# verify Amazon auth codes are available
#
[[ -z 'AWS_ACCESS_KEY_ID' ]] && { echo 'Missing $AWS_ACCESS_KEY_ID variable'; exit 1; }
[[ -z 'AWS_SECRET_ACCESS_KEY' ]] && { echo 'Missing $AWS_SECRET_KEY variable'; exit 1; }

mkdir -p $HOME/bin
cat <<EOF >> $HOME/.bash_profile
export AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY
EOF

#
# mount instance storage
#
if [[ $X_TYPE == r3.* && $X_DEVICE && $X_MOUNT ]]; then
sudo mkfs.ext4 -E nodiscard $X_DEVICE
sudo mkdir -p $X_MOUNT
sudo mount -o discard $X_DEVICE $X_MOUNT
sudo chown -R ec2-user:ec2-user $X_MOUNT
sudo chmod 775 $X_MOUNT
# set nextflow scratch path to the mounted storage
export NXF_TEMP=$X_MOUNT
echo "export NXF_TEMP=$X_MOUNT" >> $HOME/.bash_profile
fi

#
# Update docker
#
if [[ $DOCKER_VERSION ]]; then
sudo service docker stop
sudo wget https://get.docker.com/builds/Linux/x86_64/docker-${DOCKER_VERSION} -q -O /usr/bin/docker
sudo chmod +x /usr/bin/docker
sudo service docker start
fi

#
# Launch docker and pull the container when DOCKER variable is defined
#
[[ "$DOCKER_IMAGE" ]] && for x in $DOCKER_IMAGE; do docker pull $x || true; done

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

#
# Download and install nextflow
#
curl -fsSL http://www.nextflow.io/releases/${version}/nextflow  > $HOME/nextflow
chmod +x $HOME/nextflow

# launch the nextflow daemon
if [[ $NXF_ROLE != master ]]; then 
  bash -x $HOME/nextflow node -bg -cluster.join "s3:$AWS_S3BUCKET" -cluster.interface eth0
fi 

# save the environment for debugging 
env | sort > boot.env

# pull the nextflow pipeline repo
[[ $NXF_PULL ]] && $HOME/nextflow pull "$NXF_PULL" 
