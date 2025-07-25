(aws-page)=

# Amazon Web Services

:::{tip}
This page describes how to manually set up and use Nextflow with AWS Cloud.
You may be interested in using [Batch Forge](https://docs.seqera.io/platform/latest/compute-envs/aws-batch) in [Seqera Platform](https://seqera.io/platform/),
which automatically creates the required AWS infrastructure for you with minimal intervention.
:::

## AWS security credentials

Nextflow uses the [AWS security credentials](https://docs.aws.amazon.com/general/latest/gr/aws-sec-cred-types.html) to make programmatic calls to AWS services.

The AWS credentials are selected from the following sources, in order of descending priority:

1. Nextflow configuration file - `aws.accessKey` and `aws.secretKey`. See {ref}`AWS configuration<config-aws>` for more details.

2. A custom profile in `$HOME/.aws/credentials` and/or `$HOME/.aws/config`. The profile can be supplied from the `aws.profile` config option, or the `AWS_PROFILE` or `AWS_DEFAULT_PROFILE` environmental variables.

3. Environment variables - `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`.

4. The `default` profile in `~/.aws/credentials` and/or `~/.aws/config`.

5. Single Sign-On (SSO) credentials. See the [AWS documentation](https://docs.aws.amazon.com/cli/latest/userguide/sso-configure-profile-token.html) for more details.

   :::{versionadded} 23.07.0-edge
   :::

6. EC2 instance profile credentials. See the [AWS documentation](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html) and [this blog post](https://aws.amazon.com/blogs/security/granting-permission-to-launch-ec2-instances-with-iam-roles-passrole-permission/) for more details.

The AWS region is selected from the following sources, in order of descending priority:

1. Nextflow configuration file - `aws.region`
2. Environment variables - `AWS_REGION` or `AWS_DEFAULT_REGION`
3. EC2 instance metadata (if Nextflow is running in an EC2 instance)

SSO credentials and instance profile credentials are the most recommended because they don't require you to manage and distribute AWS keys explicitly. SSO credentials are ideal for launching pipelines from outside of AWS (e.g. your laptop), while instance profile credentials are ideal for launching pipelines within AWS (e.g. an EC2 instance).

## AWS IAM policies

[IAM policies](https://docs.aws.amazon.com/IAM/latest/UserGuide/access_policies.html) are the mechanism used by AWS to defines permissions for IAM identities. In order to access certain AWS services, the proper policies must be attached to the identity associated to the AWS credentials.

Minimal permissions policies to be attached to the AWS account used by Nextflow are:

- To use AWS Batch:

  ```json
  "batch:CancelJob"
  "batch:DescribeComputeEnvironments"
  "batch:DescribeJobDefinitions"
  "batch:DescribeJobQueues"
  "batch:DescribeJobs"
  "batch:ListJobs"
  "batch:RegisterJobDefinition"
  "batch:SubmitJob"
  "batch:TagResource"
  "batch:TerminateJob"
  ```

- To view [EC2](https://aws.amazon.com/ec2/) instances:

  ```json
  "ec2:DescribeInstanceAttribute"
  "ec2:DescribeInstances"
  "ec2:DescribeInstanceStatus"
  "ec2:DescribeInstanceTypes"
  "ecs:DescribeContainerInstances"
  "ecs:DescribeTasks"
  ```

- To pull container images from [ECR](https://aws.amazon.com/ecr/) repositories:

  ```json
  "ecr:BatchCheckLayerAvailability"
  "ecr:BatchGetImage"
  "ecr:DescribeImages"
  "ecr:DescribeImageScanFindings"
  "ecr:DescribeRepositories"
  "ecr:GetAuthorizationToken"
  "ecr:GetDownloadUrlForLayer"
  "ecr:GetLifecyclePolicy"
  "ecr:GetLifecyclePolicyPreview"
  "ecr:GetRepositoryPolicy"
  "ecr:ListImages"
  "ecr:ListTagsForResource"
  ```

  Alternatively, you can use AWS provided `AmazonEC2ContainerRegistryReadOnly` managed policy.

:::{note}
If you are running Fargate or Fargate Spot, you may need the following policies in addition to the listed above:
  ```json
  "ec2:DescribeSubnets"
  "ecs:CreateCluster"
  "ecs:DeleteCluster"
  "ecs:DescribeClusters"
  "ecs:ListClusters"
  ```
:::

### S3 policies

Nextflow also requires policies to access [S3 buckets](https://aws.amazon.com/s3/) in order to use the work directory, pull input data, and publish results.

Depending on the pipeline configuration, the above actions can be done all in a single bucket but, more likely, spread across multiple buckets. Once the list of buckets used by the pipeline is identified, there are two alternative ways to give Nextflow access to these buckets:

1. Grant access to all buckets by attaching the policy `"s3:*"` to the IAM identity. This works only if buckets do not set their own access policies (see point 2);

2. For more fine grained control, assign to each bucket the following policy (replace the placeholders with the actual values):

    ```json
    {
        "Version": "2012-10-17",
        "Id": "<my policy id>",
        "Statement": [
            {
                "Sid": "<my statement id>",
                "Effect": "Allow",
                "Principal": {
                    "AWS": "<ARN of the nextflow identity>"
                },
                "Action": [
                    "s3:GetObject",
                    "s3:PutObject",
                    "s3:DeleteObject",
                    "s3:PutObjectTagging",
                    "s3:AbortMultipartUpload"
                ],
                "Resource": "arn:aws:s3:::<bucket name>/*"
            },
            {
                "Sid": "AllowSSLRequestsOnly",
                "Effect": "Deny",
                "Principal": "*",
                "Action": "s3:*",
                "Resource": [
                    "arn:aws:s3:::<bucket name>",
                    "arn:aws:s3:::<bucket name>/*"
                ],
                "Condition": {
                    "Bool": {
                        "aws:SecureTransport": "false"
                    }
                }
            }
        ]
    }
    ```

See the [bucket policy documentation](https://docs.aws.amazon.com/config/latest/developerguide/s3-bucket-policy.html) for additional details.

(aws-batch)=

## AWS Batch

[AWS Batch](https://aws.amazon.com/batch/) is a managed computing service that allows the execution of containerised workloads in the AWS cloud infrastructure. It dynamically provisions the optimal quantity and type of compute resources (e.g., CPU or memory optimized compute resources) based on the volume and specific resource requirements of the jobs submitted.

Nextflow provides built-in support for AWS Batch, allowing the seamless deployment of Nextflow pipelines in the cloud, in which tasks are offloaded as Batch jobs.

Read the {ref}`AWS Batch executor <awsbatch-executor>` section to learn more about the `awsbatch` executor in Nextflow.

(aws-batch-cli)=

### AWS CLI

:::{tip}
The need for the AWS CLI is considered a legacy requirement for the deployment of Nextflow pipelines with AWS Batch.
Instead, consider using {ref}`wave-page` and {ref}`fusion-page` to facilitate access to S3 without using the AWS CLI.
:::

Nextflow uses the [AWS command line tool](https://aws.amazon.com/cli/) (`aws`) to stage input files and output files between S3 and the task containers.

The `aws` command can be made available by either (1) installing it in the container image(s) or (2) installing it in a {ref}`custom AMI <aws-custom-ami>` to be used instead of the default AMI when configuring AWS Batch.

### Get started

1. In the AWS Console, navigate to **AWS Batch** and create a [Compute environment](http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html) (CE).

   1. If you are using a custom AMI (see following sections), the AMI ID must be specified in the CE configuration
   2. Make sure to select an AMI (either custom or existing) with Docker installed (see following sections)
   3. Make sure the policy `AmazonS3FullAccess` (granting access to S3 buckets) is attached to the instance role configured for the CE
   4. If you plan to use Docker images from Amazon ECS container, make sure the `AmazonEC2ContainerServiceforEC2Role` policy is also attached to the instance role

2. In the AWS Console, create (at least) one [Job Queue](https://docs.aws.amazon.com/batch/latest/userguide/job_queues.html) and bind it to the Compute environment.

3. In the AWS Console, create an S3 bucket for the work directory (see below). You can also create separate buckets for input data and results, as needed.

4. Make sure that every process in your pipeline specifies a Docker container with the {ref}`process-container` directive.

5. Make sure that all of your container images are published in a Docker registry that can be reached by AWS Batch, such as [Docker Hub](https://hub.docker.com/), [Quay](https://quay.io/), or [Elastic Container Registry](https://aws.amazon.com/ecr/).

### Configuration

To configure your pipeline for AWS Batch:

1. Specify the AWS Batch {ref}`executor <awsbatch-executor>`
2. Specify the AWS Batch queue with the {ref}`process-queue` directive
3. Specify any Batch job container options with the {ref}`process-containerOptions` directive.

An example `nextflow.config` file is shown below:

```groovy
process {
    executor = 'awsbatch'
    queue = 'my-batch-queue'
    container = 'quay.io/biocontainers/salmon'
    containerOptions = '--shm-size 16000000 --ulimit nofile=1280:2560 --ulimit nproc=16:32'
}

aws {
    batch {
        // NOTE: this setting is only required if the AWS CLI is installed in a custom AMI
        cliPath = '/home/ec2-user/miniconda/bin/aws'
    }
    region = 'us-east-1'
}
```

:::{tip}
Each process can be configured with its own queue by using the {ref}`process-queue` directive in the process definition or via {ref}`config-process-selectors` in your Nextflow configuration.
:::

## Container Options

:::{versionadded} 21.12.1-edge
:::

The {ref}`process-containerOptions` directive can be used to control the properties of the container execution associated with each Batch job.

The following container options are currently supported:

```
-e, --env string
    Set environment variables (format: <name> or <name>=<value>)
--init
    Run an init inside the container that forwards signals and reaps processes
--memory-swap int
    The total amount of swap memory (in MiB) the container can use: '-1' to enable unlimited swap
--memory-swappiness int
    Tune container memory swappiness (0 to 100) (default -1)
--privileged
    Give extended privileges to the container
--read-only
    Mount the container's root filesystem as read only
--shm-size int
    Size (in MiB) of /dev/shm
--tmpfs string
    Mount a tmpfs directory (format: <path>:<options>,size=<int>), size is in MiB
-u, --user string
    Username or UID (format: <name|uid>[:<group|gid>])
--ulimit string
    Ulimit options (format: <type>=<soft limit>[:<hard limit>])
```

Container options may be passed in long form (e.g `--option value`) or short form (e.g. `-o value`) where available.

Few examples:

```nextflow
containerOptions '--tmpfs /run:rw,noexec,nosuid,size=128 --tmpfs /app:ro,size=64'

containerOptions '-e MYVAR1 --env MYVAR2=foo2 --env MYVAR3=foo3 --memory-swap 3240000 --memory-swappiness 20 --shm-size 16000000'

containerOptions '--ulimit nofile=1280:2560 --ulimit nproc=16:32 --privileged'
```

Check the [AWS documentation](https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html) for further details.

(aws-custom-ami)=

## Custom AMI

There are several reasons why you might need to create your own [AMI (Amazon Machine Image)](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html) to use in your Compute Environments:

- You do not want to install the AWS CLI into each of your Docker images and would rather provide it through the AMI
- The existing AMI (selected from the marketplace) does not have Docker installed
- You need to attach more storage to your EC2 instance (the default ECS instance AMI has only a 30GB EBS volume which is not enough for most data pipelines)
- You need to install additional software that is not available in your Docker image

### Create your custom AMI

From the EC2 Dashboard, select **Launch Instance**, then select **Browse more AMIs**. In the new page, select
**AWS Marketplace AMIs**, and then search for `Amazon ECS-Optimized Amazon Linux 2 (AL2) x86_64 AMI`. Select the AMI and continue as usual to configure and launch the instance.

:::{note}
The selected instance has a root volume of 30GB. Make sure to increase its size or add a second EBS volume with enough storage for real genomic workloads.
:::

When the instance is running, SSH into it (or connect with the Session Manager service), install the AWS CLI, and install any other tool that may be required (see following sections).

Finally, select **Create Image** from the EC2 Dashboard to create a new AMI from the running instance (you can also do it through the AWS CLI).

The new AMI ID needs to be specified when creating the Batch Compute Environment.

:::{warning}
Any additional software must be installed on the EC2 instance *before* creating the AMI.
:::

(id2)=

### AWS CLI installation

:::{tip}
The need for the AWS CLI is considered a legacy requirement for the deployment of Nextflow pipelines with AWS Batch.
Instead, consider using {ref}`wave-page` and {ref}`fusion-page` to facilitate access to S3 without using the AWS CLI.
:::

The [AWS CLI](https://aws.amazon.com/cli) should be installed in your custom AMI using a self-contained package manager such as [Conda](https://conda.io). That way, you can control which version of Python is used by the AWS CLI (which is written in Python).

If you don't use Conda, the `aws` command will attempt to use the version of Python that is installed in the container, and it won't be able to find the necessary dependencies.

The following snippet shows how to install AWS CLI with [Miniconda](https://conda.io/miniconda.html) in the home folder:

```bash
cd $HOME
sudo yum install -y bzip2 wget
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
$HOME/miniconda/bin/conda install -c conda-forge -y awscli
rm Miniconda3-latest-Linux-x86_64.sh
```

Afterwards, verify that the AWS CLI package works correctly:

```console
$ ./miniconda/bin/aws --version
aws-cli/1.29.20 Python/3.11.4 Linux/4.14.318-241.531.amzn2.x86_64 botocore/1.31.20
```

:::{note}
The `aws` tool will be placed in a directory named `bin` in the main installation folder. Modifying this directory structure after the tool is installed will cause it to not work properly.
:::

To configure Nextflow to use this installation, specify the `aws.batch.cliPath` option in the Nextflow configuration as shown below:

```groovy
aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
```

Replace the path above with the one matching the location where the `aws` tool is installed in your AMI.

:::{versionchanged} 19.07.0
The `executor.awscli` config option was replaced by `aws.batch.cliPath`.
:::

:::{warning}
The grandparent directory of the `aws` tool will be mounted into the container at the same path as the host, e.g. `/home/ec2-user/miniconda`, which will shadow existing files in the container. Make sure you use a path that is not already present in the container.
:::

### Docker installation

Docker is required by Nextflow to execute tasks on AWS Batch. The **Amazon ECS-optimized Amazon Linux 2023 AMI** has Docker installed, however, if you create your AMI from a different AMI that does not have Docker installed, you will need to install it manually.

The following snippet shows how to install Docker on an Amazon EC2 instance:

```bash
# install Docker
sudo yum update -y
sudo yum install docker

# start the Docker service
sudo service docker start

# empower your user to run Docker without sudo
sudo usermod -a -G docker ec2-user
```

You must logging out and logging back in again to use the new `ec2-user` permissions.

These steps must be done *before* creating the AMI from the current EC2 instance.

### Amazon ECS container agent installation

The [ECS container agent](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ECS_agent.html) is a component of Amazon Elastic Container Service (Amazon ECS) and is responsible for managing containers on behalf of ECS. AWS Batch uses ECS to execute containerized jobs, therefore it requires the agent to be installed on EC2 instances within your Compute Environments.

The ECS agent is included in the **Amazon ECS-optimized Amazon Linux 2023 AMI** . If you use a different base AMI, you can also install the agent on any EC2 instance that supports the Amazon ECS specification.

To install the agent, follow these steps:

```bash
sudo yum install ecs-init
sudo systemctl enable --now ecs
```

To test the installation:

```bash
curl -s http://localhost:51678/v1/metadata | python -mjson.tool (test)
```

:::{note}
The `AmazonEC2ContainerServiceforEC2Role` policy must be attached to the instance role in order to be able to connect the EC2 instance created by the Compute Environment to the ECS container.
:::

:::{note}
The `AmazonEC2ContainerRegistryReadOnly` policy should be attached to the instance role in order to get read-only access to Amazon EC2 Container Registry repositories.
:::

## Jobs & Execution

### Custom job definition

Nextflow automatically creates the Batch [Job definitions](http://docs.aws.amazon.com/batch/latest/userguide/job_definitions.html) needed to execute tasks in your pipeline, so you don't need to define them beforehand.

However, sometimes you may still need to specify a custom **Job Definition** to fine tune the configuration of a specific job, for example to define custom mount paths.

To do that, first create a **Job Definition** in the AWS Console (or by other means). Note the name of the Job definition you created. You can then associate a process execution with this Job definition by using the {ref}`process-container` directive and specifying, in place of the container image name, the Job definition name prefixed by `job-definition://`, as shown below:

```groovy
process.container = 'job-definition://your-job-definition-name'
```

### Hybrid workloads

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment of hybrid workloads in which some jobs are executed in the local computer or local computing cluster and some jobs are offloaded to AWS Batch.

To enable this feature, use one or more {ref}`config-process-selectors` in your Nextflow configuration to apply the AWS Batch configuration to the subset of processes that you want to offload. For example:

```groovy
process {
    withLabel: bigTask {
        executor = 'awsbatch'
        queue = 'my-batch-queue'
        container = 'my/image:tag'
    }
}

aws {
    region = 'eu-west-1'
}
```

With the above configuration, processes with the `bigTask` {ref}`process-label` will run on AWS Batch, while the remaining processes will run in the local computer.

Then launch the pipeline with the -bucket-dir option to specify an AWS S3 path for the jobs computed with AWS Batch and, optionally, the -work-dir to specify the local storage for the jobs computed locally:

```bash
nextflow run <script or project name> -bucket-dir s3://my-bucket/some/path
```

:::{warning}
The AWS S3 path needs to contain at least one sub-directory (e.g. `s3://my-bucket/work` rather than `s3://my-bucket`).
:::

:::{note}
Nextflow will automatically manage the transfer of input and output files between the local and cloud environments when using hybrid workloads.
:::

### Volume mounts

:::{versionadded} 19.07.0
:::

User provided container volume mounts can be provided as shown below:

```groovy
aws {
    region = 'eu-west-1'
    batch {
        volumes = '/tmp'
    }
}
```

Multiple volumes can be specified as a comma-separated list of paths. The usual Docker volume mount syntax can be used to specify complex volumes where the container path is different from the host path or the volume should be *read-only*. For example:

```groovy
aws {
    region = 'eu-west-1'
    batch {
        volumes = ['/tmp', '/host/path:/mnt/path:ro']
    }
}
```

The above snippet defines two volume mounts for the jobs executed in your pipeline. The first volume mounts the host path `/tmp` to the same path in the container, with the *read-write* access mode. The second volume mounts the host path `/host/path` to `/mnt/path` in the container, with the *read-only* access mode.

### Troubleshooting

**Problem**: The Pipeline execution terminates with an AWS error message similar to the one shown below:

```
JobQueue <your queue> not found
```

Make sure you have defined a AWS region in the Nextflow configuration file and it matches the region in which your Batch environment has been created.

**Problem**: A process execution fails reporting the following error message:

```
Process <your task> terminated for an unknown reason -- Likely it has been terminated by the external system
```

This may happen when Batch is unable to execute the process script. A common cause of this problem is that the Docker container image you have specified uses a non standard [entrypoint](https://docs.docker.com/engine/reference/builder/#entrypoint) which does not allow the execution of the Bash launcher script required by Nextflow to run the job.

This may also happen if the AWS CLI doesn't run correctly.

Other places to check for error information:

- The `.nextflow.log` file.
- The Job execution log in the AWS Batch dashboard.
- The [CloudWatch](https://aws.amazon.com/cloudwatch/) logs found in the `/aws/batch/job` log group.

**Problem**: A process execution is stalled in the `RUNNABLE` status and the pipeline output is similar to the one below:

```
executor >  awsbatch (1)
process > <your process> (1) [  0%] 0 of ....
```

It may happen that the pipeline execution hangs indefinitely because one of the jobs is held in the queue and never gets executed. In AWS Console, the queue reports the job as `RUNNABLE` but it never moves from there.

There are multiple reasons why this can happen. They are mainly related to the Compute Environment workload/configuration, the docker service or container configuration, network status, etc.

This [AWS page](https://aws.amazon.com/premiumsupport/knowledge-center/batch-job-stuck-runnable-status/) provides several resolutions and tips to investigate and work around the issue.

(aws-fargate)=

## AWS Fargate

:::{versionadded} 23.12.0-edge
:::

Nextflow provides experimental support for the execution of [AWS Batch jobs with Fargate resources](https://docs.aws.amazon.com/batch/latest/userguide/fargate.html).

AWS Fargate is a technology that you can use with AWS Batch to run containers without having to manage servers or  EC2 instances.
With AWS Fargate, you no longer have to provision, configure, or scale clusters of virtual machines to run containers.

To enable the use of AWS Fargate in your pipeline use the following settings in your `nextflow.config` file:

```groovy
process.executor = 'awsbatch'
process.queue = '<AWS BATCH QUEUE>'
aws.region = '<AWS REGION>'
aws.batch.platformType = 'fargate'
aws.batch.jobRole = 'JOB ROLE ARN'
aws.batch.executionRole = 'EXECUTION ROLE ARN'
wave.enabled = true
```

See the AWS documentation for details how to create the required AWS Batch queue for Fargate, the Batch Job Role
and the Batch Execution Role.

:::{note}
Nextflow uses [s5cmd](https://github.com/peak/s5cmd) to download the task input data and upload the task outputs.
To enable this capability, you need to enable the Wave service in the Nextflow configuration, as shown in the above example. See {ref}`Wave <wave-page>` documentation for more details.
:::

## Advanced configuration

Read the {ref}`AWS configuration<config-aws>` section to learn more about advanced configuration options.
