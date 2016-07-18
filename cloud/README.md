# Cloud boot 

* The script `launch.sh` launches a set of EC2 nodes and the master node. 
    
* Each EC2 node on start downloads and runs the script `cloud-boot.sh`. This script does the following: 
    - Downloads the Nextflow installer script
    - Pulls the Docker image specified by the `DOCKER_IMAGE` environment variable (multiple images by using a blank character separated list).
    - Launches the Nextflow daemon in background
    - Tag EC2 instances with a `master`/`worker` label
    - Mount instance ephemeral storage 
    - Mount [EFS](https://aws.amazon.com/efs/) storage 
    
* It requires a Linux Amazon image with the following tools installed:
    - cloud-init
    - curl 
    - docker 
    - java se 7 or 8
    
NOTE: if you use an AMI other than Amazon Linux, change the user name `ec2-user` in the script `cloud-boot.sh` 
accordingly. 
    
    
## Quickstart
    
* Define your configuration details setting properly the variables in the `launch.sh` script, in particular: 
    
    - AWS_ACCESS_KEY_ID=...
    - AWS_SECRET_ACCESS_KEY=...
    - X_AMI=<image ID eg. ami-f71a7f80>
    - X_TYPE=<instance type e.g. m3.xlarge>
    - X_SECURITY=<security group e.g sg-72b74a05>
    - X_KEY=<security key name>
    - X_BUCKET=<s3 bucket where nodes share their IPs in order to discover each other>
  
    
* Use the script `launch.sh` to setup the cluster, for example: 
 
        $ ./launch cluster 10 
       
* Add new worker nodes with the `node` command. For example:

        $ ./launch.sh node 5 
      
* SSH into the master node

* Launch at nextflow pipeline execution     
