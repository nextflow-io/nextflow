#!/bin/sh
set -e
set -u

#
# DEFINE the following variables
#

#X_RUN=<unique id of the execution>

#AWS_ACCESS_KEY_ID=xxx
#AWS_SECRET_ACCESS_KEY=yyy
#X_AMI=<image ID eg. ami-f71a7f80>
#X_TYPE=<instance type e.g. m3.xlarge>
#X_SECURITY=<security group e.g sg-72b74a05>
#X_KEY=<security key name>
#X_BUCKET=<s3 bucket where nodes share their IPs in order to discover each other>
#X_SUBNET=<vpc subnet id>

# only for spot instances
#X_ZONE=eu-west-1c
#X_PRICE=0.125

# boot storage size
#X_STORAGE

# only required to mount instance storage
#X_DEVICE=/dev/xvdc
#X_MOUNT=/mnt/scratch

# EFS file system
#X_EFS_ID=fs-xxx
#X_EFS_MOUNT=/mnt/efs

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
echo "Run           : $X_RUN"
echo "AMI           : $X_AMI"
echo "Instance type : $X_TYPE"
echo "Security group: $X_SECURITY"
echo "Security key  : $X_KEY"
echo "S3 join bucket: $X_BUCKET"

set +u
echo "Root storage  : $X_STORAGE GB"
echo "EFS ID        : ${X_EFS_ID:--}"
echo "EFS mount     : ${X_EFS_MOUNT:--} "

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

#
# create the cloud boothook file
#
local boothook=$(mktemp)

if [[ $X_EFS_ID && $X_EFS_MOUNT ]]; then
cat <<EndOfString >>$boothook
zone="\$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone)"
region="\${zone::-1}"
yum install -y nfs-utils
mkdir -p $X_EFS_MOUNT
mount -t nfs4 -o nfsvers=4.1 \${zone}.${X_EFS_ID}.efs.\${region}.amazonaws.com:/ $X_EFS_MOUNT
chown ec2-user:ec2-user $X_EFS_MOUNT
chmod 775 $X_EFS_MOUNT
EndOfString
fi

if [[ $X_TYPE == r3.* && $X_DEVICE && $X_MOUNT ]]; then
cat <<EndOfString >>$boothook
mkfs.ext4 -E nodiscard $X_DEVICE
mkdir -p $X_MOUNT
mount -o discard $X_DEVICE $X_MOUNT
chown -R ec2-user:ec2-user $X_MOUNT
chmod 775 $X_MOUNT
EndOfString
fi

#
# create the script-shell  file
#
local scriptschell=$(mktemp)
cat <<EndOfString >> $scriptschell
#!/bin/bash
su - ec2-user << 'EndOfScript'
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
export X_EFS_ID="$X_EFS_ID"
export X_EFS_MOUNT="$X_EFS_MOUNT"
export X_RUN="$X_RUN"
(
$(cat ./cloud-boot.sh)
) &> ~ec2-user/boot.log
EndOfScript
EndOfString

#
# combine the two file as mime-multipart content user-data
#
if [[ -s $boothook ]]; then
  ./mimetext.py $scriptschell:x-shellscript $boothook:cloud-boothook
else
  ./mimetext.py $scriptschell:x-shellscript
fi

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

    if [[ $X_SUBNET ]]; then
    cli+=(--subnet-id); cli+=("$X_SUBNET")
    fi

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

if [[ $X_SUBNET ]]; then
local subnet="\"SubnetId\": \"$X_SUBNET\", "
fi

TMPFILE=$(mktemp)
cat >$TMPFILE <<EOF 
{
  "ImageId": "$X_AMI",
  "KeyName": "$X_KEY",
  "SecurityGroupIds": [ "$X_SECURITY" ], $subnet
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




