#!/bin/bash
# TODO: This is a temporary workaround to run nextflow on GCE cluster from development code
# After cluster has been created, nextflow will fail to start with GCE driver because it is pulling official nextflow.
# Run this after cluster has been created and the local code will be copied to the cluster and nextflow will be started on
# all worker nodes.

CLUSTER=$1
if [ -z "$CLUSTER" ]; then
  echo "Usage: $0 <cluster name>"
  ./nextflow cloud list -driver gce
  exit 1
fi
SSH_OPTS="-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -i ~/.ssh/id_rsa" 
cd $HOME
tar -cf /tmp/dot_nextflow.tar.gz .nextflow
cd -
./nextflow  cloud list $CLUSTER -driver gce | grep RUNNING | while read id add state role; do
  echo "Initializing $role $add"
  scp $SSH_OPTS /tmp/dot_nextflow.tar.gz $USER@$add:
  ssh $SSH_OPTS $USER@$add 'tar --warning=no-unknown-keyword -xf dot_nextflow.tar.gz' </dev/null
  if [[ "$role" == worker ]]; then
    echo "Starting worker process"
    ssh $SSH_OPTS $USER@$add '. ~/.bash_profile; nohup ./nextflow node -cluster.join "$NXF_CLUSTER_JOIN" -cluster.interface eth0 -bg >/dev/null 2>&1 &' </dev/null
  fi

  echo "Done"
done
