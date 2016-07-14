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

[[ $NXF_VER ]] && echo "export NXF_VER='$NXF_VER'" >> $HOME/.bash_profile
[[ $NXF_OPTS ]] && echo "export NXF_OPTS='$NXF_OPTS'" >> $HOME/.bash_profile
[[ $NXF_TRACE ]] && echo "export NXF_TRACE='$NXF_TRACE'" >> $HOME/.bash_profile

#
# set instance name
#
instance="$(curl -s http://169.254.169.254/latest/meta-data/instance-id)"
zone="$(curl -s 169.254.169.254/latest/meta-data/placement/availability-zone)"
region="${zone::-1}"
aws ec2 --region "$region" create-tags --resources "$instance" --tags "Key=Name,Value=$NXF_ROLE"
aws ec2 --region "$region" create-tags --resources "$instance" --tags "Key=Cluster,Value=$X_RUN"

#
# mount instance storage
#
if [[ $X_DEVICE && $X_MOUNT ]]; then
# set nextflow scratch path to the mounted storage
export NXF_TEMP=$X_MOUNT
echo "export NXF_TEMP=$X_MOUNT" >> $HOME/.bash_profile
fi

#
# Update docker (see https://get.docker.com/builds/ for details)
#
if [[ $DOCKER_VERSION ]]; then
sudo service docker stop
sudo wget https://get.docker.com/builds/Linux/x86_64/docker-${DOCKER_VERSION} -q -O /usr/bin/docker
sudo chmod +x /usr/bin/docker
# see https://github.com/docker/docker/issues/18113
sudo rm -rf /var/lib/docker/network/files/local-kv.db
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

# pull the nextflow pipeline repo
[[ $NXF_PULL ]] && $HOME/nextflow pull "$NXF_PULL"

# launch the nextflow daemon
cluster_join="s3:$AWS_S3BUCKET"
cluster_work="s3://<your-bucket/work>"
if [[ $X_EFS_ID && $X_EFS_MOUNT ]]; then
  [[ $X_RUN ]] && cluster_join="path:$X_EFS_MOUNT/cluster/$X_RUN"
  cluster_work=$X_EFS_MOUNT/work
fi

if [[ $NXF_ROLE == worker ]]; then
  bash -x $HOME/nextflow node -bg -cluster.join "$cluster_join" -cluster.interface eth0
else
cat <<EOF >> $HOME/README
#
# Launch the pipeline execution by using the following command
#
./nextflow run cbcrg/kallisto-nf \
  -process.executor ignite \
  -cluster.join $cluster_join \
  -cluster.interface=eth0  \
  -work-dir ${cluster_work} \
  -with-docker
EOF
fi

# save the environment for debugging 
env | sort > boot.env

# just a marker file
touch .READY
