Cloud boot 
============

* The script `launch.sh` launch a set of EC2 nodes. 

* The above script uses the `ec2-userdata.txt` to set the environment in the started nodes, including: 
    - AWS_ACCESS_KEY
    - AWS_SECRET_KEY
    - NXF_VER
    - DOCKER_IMAGE
    
* Each EC2 node on bootstrap will download and run the script `cloud-boot.sh`. This script do the following: 
    - Download the Nextflow installer script
    - Pull the Docker image specified by the `DOCKER_IMAGE` environment variable
    - Launch the Nextflow daemon in background
    
    
    
Quickstart
------------
    
* Use the script `launch.sh` to launch a set of workers 

* Check the node availability listing the content of the bucket `nxf-cluster` with the command: 
    
    aws s3 ls s3://nxf-cluster
       
* Launch the master node (modify the launcher script) 
      
* SSH into the master node

* Define the these environment variables in the master 

    - AWS_ACCESS_KEY=xxx
    - AWS_SECRET_KEY=xxx
    - NXF_VER=xx
    
* Download nextflow with the command: 

    curl -fsSL http://www.nextflow.io/releases/v${NXF_VER}/nextflow | bash