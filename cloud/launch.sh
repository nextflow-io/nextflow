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

# boot storage size
#X_STORAGE

# only required to mount instance storage
#X_DEVICE=/dev/xvdc
#X_MOUNT=/mnt/scratch

# nextflow details
#NXF_VER=0.17.3
#NXF_MODE=ignite
#NXF_OPTS='-Xms512M -Xmx2G'
#NXF_PULL=<github repo to pull>

# docker containers to pull on start
#DOCKER_IMAGE="busybox"

# docker runtime to be used
#DOCKER_VERSION=1.10.3

# Print the current status
echo "AMI           : $X_AMI"
echo "Instance type : $X_TYPE"
echo "Security group: $X_SECURITY"
echo "Security key  : $X_KEY"
echo "S3 join bucket: $X_BUCKET"

set +u
echo "Root storage  : $X_STORAGE GB"

if [[ $3 == '--spot' ]]; then
  X_SPOT=true 
  X_PRICE=${4:-$X_PRICE}
  [[ ! $X_PRICE ]] && echo "ERROR -- Spot instance bid price need to be specified" && exit
echo "Instance type : **SPOT**"
echo "Avail zone    : $X_ZONE"
echo "Spot price    : $X_PRICE"
else
echo "Instance type : **ON-DEMAND**"
fi

if [[ $1 == master ]]; then
    X_MSG="This is going to launch the MASTER node. Is it OK (y/n)?"
    X_COUNT=1
elif [[ $1 == cluster ]]; then
    X_COUNT=${2:-1}
    X_MSG="This is going to launch * $X_COUNT * instances. Is it OK (y/n)?"
elif [[ $1 == node ]]; then
    X_COUNT=${2:-1}
    X_MSG="This is going to launch * $X_COUNT * instances. Is it OK (y/n)?"
else
    echo "ERROR -- Unknown command: $1" && exit 1
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
export AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY"
export AWS_S3BUCKET="$AWS_S3BUCKET"
export DOCKER_IMAGE="$DOCKER_IMAGE"
export DOCKER_VERSION="${DOCKER_VERSION}"
export NXF_VER="$NXF_VER"
export NXF_MODE=ignite
export NXF_OPTS="$NXF_OPTS"
export NXF_PULL="$NXF_PULL"
export NXF_ROLE="$role"
export NXF_TRACE="$NXF_TRACE"
export X_TYPE="$X_TYPE"
export X_MOUNT="$X_MOUNT"
export X_DEVICE="$X_DEVICE"
curl -fsSL https://raw.githubusercontent.com/nextflow-io/nextflow/master/cloud/cloud-boot.sh | bash &> ~ec2-user/boot.log
EOF
EndOfString

}

function getRootDevice() {
  local ami=$1
  local size=$2

  local str="$(aws ec2 describe-images --image-ids $ami --query 'Images[*].{ID:BlockDeviceMappings}' --output text)"
  local device=$(echo "$str" | grep ID | cut -f 2)
  local delete=$(echo "$str" | grep EBS | cut -f 2 | tr '[:upper:]' '[:lower:]')
  local snapsh=$(echo "$str" | grep EBS | cut -f 4)
  local type=$(echo "$str" | grep EBS | cut -f 6)

cat << EndOfString
{
    "DeviceName": "$device",
    "Ebs": {
        "DeleteOnTermination": $delete,
        "SnapshotId": "$snapsh",
        "VolumeSize": $size,
        "VolumeType": "$type"
    }
}
EndOfString

}

function getDeviceMapping() {

if [[ $X_DEVICE ]]; then
cat << EndOfString
{
    "DeviceName": "$X_DEVICE",
    "VirtualName": "ephemeral0"
}
EndOfString
fi

if [[ $X_STORAGE ]]; then
[[ $X_DEVICE ]] && echo ","
getRootDevice $X_AMI $X_STORAGE
fi

}

#
# launch EC2 on-request instances
#
function runInstances() {

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

    if [[ $X_DEVICE || $X_STORAGE ]]; then
    cli+=(--block-device-mappings); cli+=("[$(getDeviceMapping)]")
    fi

    "${cli[@]}"

}

#
# Launches spot instances
# 
function runSpot() {

local role=$1
local count=$2

if [[ $X_DEVICE || $X_STORAGE ]]; then
local mapping="\"BlockDeviceMappings\": [$(getDeviceMapping)], "
fi

TMPFILE=$(mktemp)
cat >$TMPFILE <<EOF 
{
  "ImageId": "$X_AMI",
  "KeyName": "$X_KEY",
  "SecurityGroupIds": [ "$X_SECURITY" ],
  "InstanceType": "$X_TYPE",
  "UserData": "$(cloudInit $role | base64)", $mapping
  "Placement": {
    "AvailabilityZone": "$X_ZONE"
  }
}
EOF
    
aws ec2 request-spot-instances \
  --spot-price "$X_PRICE" \
  --instance-count "$count" \
  --launch-specification "file://$TMPFILE"  

}


#
# launch a cluster i.e. master + nodes
#
function launch_cluster() {
  local x=$(($X_COUNT-1))

  if [[ $X_SPOT ]]; then
    runSpot master 1
    (( $x > 0 )) && runSpot worker $x
  else
    runInstances master 1
    (( $x > 0 )) && runInstances worker $x
  fi
}


#
# Launches the master node
#
function launch_master() {
    runInstances master 1
}

#
# Launches on-request instances
#
function launch_nodes() {
  if [[ $X_SPOT ]]; then
    runSpot worker $X_COUNT
  else
    runInstances worker $X_COUNT
  fi
}

# Confirmation
read -p "$X_MSG " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then echo ABORTED; exit 1; fi

if [[ $1 == cluster ]]; then
    # clean-up S3 bucket holding IP nodes to join
    aws s3 rm "s3://$X_BUCKET/" --recursive
    launch_cluster
elif [[ $1 == master ]]; then
    launch_master
else
    launch_nodes
fi




