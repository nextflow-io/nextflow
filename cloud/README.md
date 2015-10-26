Cloud boot 
============

* The script `launch.sh` to launch a set of EC2 nodes and the master node. 

* The above script uses the `ec2-userdata.txt` to set the environment in the started nodes, including: 
    - AWS_ACCESS_KEY
    - AWS_SECRET_KEY
    - NXF_VER
    - DOCKER_IMAGE
    
* Each EC2 node on start will download and run the script `cloud-boot.sh`. This script do the following: 
    - Download the Nextflow installer script
    - Pull the Docker image specified by the `DOCKER_IMAGE` environment variable
    - Launch the Nextflow daemon in background
    
* This scripts require a Linux Amazon image having at least the following tools installed:
    - cloud-init
    - curl 
    - docker 
    
NOTE: if you use an AMI other than Amazon Linux, change the user name `ec2-user` in the script `cloud-boot.sh` 
accordingly. 
    
Quickstart
------------
    
* Use the script `launch.sh` to launch a set of workers: 
 
    $ ./launch 5 

* Check the node availability listing the content of the bucket `nxf-cluster` with the command: 
    
    $ aws s3 ls s3://nxf-cluster
       
* Launch the master node:

    $ ./launch.sh master
      
* SSH into the master node

* Install nextflow with this command

    $ curl -fsSL http://get.nextflow.io | bash

* Define the AWS environment variables in the master 

    - AWS_ACCESS_KEY=xxx
    - AWS_SECRET_KEY=xxx
    
* Launch at nextflow pipeline execution     