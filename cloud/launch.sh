#!/bin/sh
set -e
set -u

#
# DEFINE the following variables
#
#X_AMI=<image ID eg. ami-f71a7f80>
#X_TYPE=${instance type e.g. m3.xlarge}
#X_SECURITY=<security group e.g sg-72b74a05>
#X_KEY=<security key name>
#X_BUCKET=<s3 bucket where nodes share their IPs in order to discover each other>
#X_USERDATA=<user-data file e.g. ec2-userdata.txt>
#

# Print the current status
echo "AMI           : $X_AMI"
echo "Instance type : $X_TYPE"
echo "Security group: $X_SECURITY"
echo "Security key  : $X_KEY"
echo "S3 join bucket: $X_BUCKET"
echo "User-data     : $X_USERDATA"
echo ''

set +u
if [[ $1 == master ]]; then
X_MSG="This is going to launch the MASTER node. Is it OK (y/n)?"
else
X_COUNT=${1:-1}
X_MSG="This is going to launch * $X_COUNT * instance. Is it OK (y/n)?"
fi

function launch_master() {
aws ec2 run-instances \
    --image-id "$X_AMI" \
    --instance-type "$X_TYPE" \
    --security-group-ids "$X_SECURITY" \
    --key-name "$X_KEY" \
    --count 1

}

function launch_nodes() {
aws ec2 run-instances \
    --image-id "$X_AMI" \
    --instance-type "$X_TYPE" \
    --security-group-ids "$X_SECURITY" \
    --key-name "$X_KEY" \
    --count "$X_COUNT" \
    --user-data "file:$X_USERDATA"
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
launch_nodes
fi




