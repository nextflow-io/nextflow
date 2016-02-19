#!/bin/sh
set -e
set -u

#
# DEFINE the following variables
#
#AWS_ACCESS_KEY_ID=xxx
#AWS_SECRET_ACCESS_KEY=yyy
#X_AMI=<image ID eg. ami-f71a7f80>
#X_TYPE=<instance type e.g. m3.xlarge>
#X_SECURITY=<security group e.g sg-72b74a05>
#X_KEY=<security key name>
#X_BUCKET=<s3 bucket where nodes share their IPs in order to discover each other>

# only for spot instances
#X_ZONE=eu-west-1c
#X_PRICE=0.125

# only required to mount instance storage
#X_DEVICE=/dev/xvdc
#X_MOUNT=/mnt/scratch

# nextflow details
#NXF_VER=0.17.3
#NXF_MODE=ignite
#NXF_PULL=<github repo to pull>

# docker containers to pull on start
#DOCKER_IMAGE="grape/contig:rgcrg-0.1 grape/quantification:flux-1.6.1 grape/inferexp:rseqc-2.3.9 grape/mapping:gem-1.7.1"

# docker runtime to be used
#DOCKER_VERSION=1.10.1

# Print the current status
echo "AMI           : $X_AMI"
echo "Instance type : $X_TYPE"
echo "Security group: $X_SECURITY"
echo "Security key  : $X_KEY"
echo "S3 join bucket: $X_BUCKET"


set +u

if [[ $2 == '--spot' ]]; then 
  X_SPOT=true 
  X_PRICE=${3:-$X_PRICE}
  [[ ! $X_PRICE ]] && echo "ERROR -- Spot instance bid price need to be specified" && exit
echo "Avail zone    : $X_ZONE"
echo "Spot price    : $X_PRICE"
fi 

if [[ $1 == master ]]; then
X_MSG="This is going to launch the MASTER node. Is it OK (y/n)?"
else
X_COUNT=${1:-1}
X_MSG="This is going to launch * $X_COUNT * instance. Is it OK (y/n)?"
fi

function get_abs_filename() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

#
# Prints the cloud-init string 
# 
function cloudInit() {

local role=$1
cat << EndOfString
#!/bin/bash
su - ec2-user << 'EOF'
export AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY
export AWS_S3BUCKET=$AWS_S3BUCKET
export DOCKER_IMAGE=$DOCKER_IMAGE
export DOCKER_VERSION=${DOCKER_VERSION}
export NXF_VER=$NXF_VER
export NXF_MODE=ignite
export NXF_PULL=$NXF_PULL
export NXF_ROLE=$role
export X_TYPE=$X_TYPE
export X_MOUNT=$X_MOUNT
export X_DEVICE=$X_DEVICE
curl -fsSL https://raw.githubusercontent.com/nextflow-io/nextflow/master/cloud/cloud-boot.sh | bash &> ~ec2-user/boot.log
EOF
EndOfString

}

#
# launch EC2 on-request instances
#
function launchInstances() {

    local role=$1
    local count=$2

    declare -a cli=()
    cli+=(aws)
    cli+=(ec2)
    cli+=(run-instances)
    cli+=(--count); cli+=($count)
    cli+=(--image-id); cli+=("$X_AMI")
    cli+=(--instance-type); cli+=("$X_TYPE")
    cli+=(--security-group-ids); cli+=("$X_SECURITY")
    cli+=(--key-name); cli+=("$X_KEY")
    cli+=(--user-data); cli+=("$(cloudInit $role)")

    if [[ $X_DEVICE ]]; then
    local mapping="[{\"DeviceName\": \"$X_DEVICE\",\"VirtualName\": \"ephemeral0\"}]"
    cli+=(--block-device-mappings); cli+=("$mapping")
    fi

    "${cli[@]}"

}

#
# Launches the master node
# 
function launch_master() {
    launchInstances master 1
}

#
# Launches on-request instances
# 
function launch_nodes() {
    launchInstances worker $X_COUNT
}

#
# Launches spot instances
# 
function launch_spot() {

if [[ $X_DEVICE ]]; then
local mapping="\"BlockDeviceMappings\": [{\"DeviceName\": \"$X_DEVICE\",\"VirtualName\": \"ephemeral0\"}], "
fi

TMPFILE=$(mktemp)
cat >$TMPFILE <<EOF 
{
  "ImageId": "$X_AMI",
  "KeyName": "$X_KEY",
  "SecurityGroupIds": [ "$X_SECURITY" ],
  "InstanceType": "$X_TYPE",
  "UserData": "$(cloudInit worker | base64)", $mapping
  "Placement": {
    "AvailabilityZone": "$X_ZONE"
  }
}
EOF
    
aws ec2 request-spot-instances \
  --spot-price "$X_PRICE" \
  --instance-count "$X_COUNT" \
  --launch-specification "file://$TMPFILE"  

}


# Confirmation
read -p "$X_MSG " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then echo ABORTED; exit 1; fi

if [[ $1 == master ]]; then
launch_master
else
# empty the cluster bucket
aws s3 rm "s3://$X_BUCKET/" --recursive
[[ $X_SPOT ]] && launch_spot || launch_nodes
fi




