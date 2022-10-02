.. _aws-page:

************
Amazon Cloud
************

AWS security credentials
=========================

Nextflow uses the `AWS security credentials <https://docs.aws.amazon.com/general/latest/gr/aws-sec-cred-types.html>`_
to make programmatic calls to AWS services.

You can provide your AWS access keys using the standard AWS variables shown below:

    * ``AWS_ACCESS_KEY_ID``
    * ``AWS_SECRET_ACCESS_KEY``
    * ``AWS_DEFAULT_REGION``

If ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY`` are not defined in the environment, Nextflow will attempt to
retrieve credentials from your ``~/.aws/credentials`` or ``~/.aws/config`` files. The ``default`` profile can be
overridden via the environmental variable ``AWS_PROFILE`` (or ``AWS_DEFAULT_PROFILE``).

Alternatively AWS credentials can be specified in the Nextflow configuration file.

See :ref:`AWS configuration<config-aws>` for more details.

.. note:: Credentials can also be provided by using an IAM Instance Role. The benefit of this approach is that
  it spares you from managing/distributing AWS keys explicitly.
  Read the `IAM Roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_ documentation
  and `this blog post <https://aws.amazon.com/blogs/security/granting-permission-to-launch-ec2-instances-with-iam-roles-passrole-permission/>`_ for more details.


AWS IAM policies
=================

`IAM policies <https://docs.aws.amazon.com/IAM/latest/UserGuide/access_policies.html>`_ are the mechanism used by AWS to
defines permissions for IAM identities. In order to access certain AWS services, the proper policies must be
attached to the identity associated to the AWS credentials.

Minimal permissions policies to be attached to the AWS account used by Nextflow are:

- To interface AWS Batch::

  "batch:DescribeJobQueues"
  "batch:CancelJob"
  "batch:SubmitJob"
  "batch:ListJobs"
  "batch:DescribeComputeEnvironments"
  "batch:TerminateJob"
  "batch:DescribeJobs"
  "batch:RegisterJobDefinition"
  "batch:DescribeJobDefinitions"

- To be able to see the `EC2 <https://aws.amazon.com/ec2/>`_ instances::

  "ecs:DescribeTasks"
  "ec2:DescribeInstances"
  "ec2:DescribeInstanceTypes"
  "ec2:DescribeInstanceAttribute"
  "ecs:DescribeContainerInstances"
  "ec2:DescribeInstanceStatus"

- To pull container images stored in the `ECR <https://aws.amazon.com/ecr/>`_ repositories::

  "ecr:GetAuthorizationToken"
  "ecr:BatchCheckLayerAvailability"
  "ecr:GetDownloadUrlForLayer"
  "ecr:GetRepositoryPolicy"
  "ecr:DescribeRepositories"
  "ecr:ListImages"
  "ecr:DescribeImages"
  "ecr:BatchGetImage"
  "ecr:GetLifecyclePolicy"
  "ecr:GetLifecyclePolicyPreview"
  "ecr:ListTagsForResource"
  "ecr:DescribeImageScanFindings"

S3 policies
------------
Nextflow requires policies also to access `S3 buckets <https://aws.amazon.com/s3/>`_ in order to:

1. use the workdir
2. pull input data
3. publish results

Depending on the pipeline configuration, the above actions can be done all in a single bucket but, more likely, spread across multiple
buckets. Once the list of buckets used by the pipeline is identified, there are two alternative ways to give Nextflow access to these buckets:

1. grant access to all buckets by attaching the policy ``"s3:*"`` to the AIM identity. This works only if buckets do not set their own access policies (see point 2);
2. for a more fine grained control, assign to each bucket the following policy (replace the placeholders with the actual values)::

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
                "s3:DeleteObject"
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

See the `bucket policy documentation <https://docs.aws.amazon.com/config/latest/developerguide/s3-bucket-policy.html>`_
for additional details.


.. _aws-batch:

AWS Batch
=========

.. note::
    Requires Nextflow version `0.26.0` or later.

`AWS Batch <https://aws.amazon.com/batch/>`_ is a managed computing service that allows the execution of containerised
workloads in the Amazon cloud infrastructure. It dynamically provisions the optimal quantity and type of compute
resources (e.g., CPU or memory optimized compute resources) based on the volume and specific resource requirements
of the jobs submitted.

Nextflow provides a built-in support for AWS Batch which allows the seamless deployment of a Nextflow pipeline
in the cloud offloading the process executions as Batch jobs.

.. _aws-batch-config:

AWS CLI
--------

Nextflow requires to access the `AWS command line tool <https://aws.amazon.com/cli/>`_ (``aws``) from the container in
which the job runs in order to stage the required input files and to copy back the resulting output files in the S3 storage.

The ``aws`` tool can be made available in the container in two ways:

1. installed in the Docker image(s) used during the pipeline execution,

2. installed in a custom `AMI (Amazon Machine Image) <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html>`_ to use in place of the default AMI when configuring AWS Batch (see next section).

The latter approach is preferred because it allows the use of existing Docker images without the need to add
the AWS CLI tool to them.

See the sections below to learn how to create a custom AMI and install the AWS CLI tool to it.

Get started
-------------

1. In the AWS Console, create a `Compute environment <http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html>`_ (CE) in your AWS Batch Service.

    a. If you are using a custom AMI (see following sections), the AMI ID must be specified in the CE configuration
    b. Make sure to select an AMI (either custom or existing) with Docker installed (see following sections)
    c. Make sure the policy ``AmazonS3FullAccess`` (granting access to S3 buckets) is attached to the instance role configured for the CE
    d. If you plan to use Docker images from Amazon ECS container, make sure the ``AmazonEC2ContainerServiceforEC2Role`` policy is also attached to the instance role

2. In the AWS Console, create (at least) one `Job Queue <https://docs.aws.amazon.com/batch/latest/userguide/job_queues.html>`_ and bind it to the Compute environment.

3. In the AWS Console, create an S3 storage's bucket for the bucket-dir (see below) and others for the input data and results, as needed.

4. Make sure your pipeline processes specifies one or more Docker containers by using the :ref:`process-container` directive.

5. Container images need to be published in a Docker registry such as `Docker Hub <https://hub.docker.com/>`_, `Quay <https://quay.io/>`_ or `ECS Container Registry <https://aws.amazon.com/ecr/>`_ that can be reached by ECS Batch.

Configuration
-------------

When configuring your pipeline:

1. import the `nf-amazon` plugin
2. specify the AWS Batch :ref:`executor<awsbatch-executor>`
3. specify one or more AWS Batch queues for the execution by using the :ref:`process-queue` directive
4. specify the AWS job container properties by using the :ref:`process-containerOptions` directive.

An example ``nextflow.config`` file is shown below::

    plugins {
        id 'nf-amazon'
    }

    process {
        executor = 'awsbatch'
        queue = 'my-batch-queue'
        container = 'quay.io/biocontainers/salmon'
        containerOptions = '--shm-size 16000000 --ulimit nofile=1280:2560 --ulimit nproc=16:32'
    }

    aws {
        batch {
            // NOTE: this setting is only required if the AWS CLI tool is installed in a custom AMI
            cliPath = '/home/ec2-user/miniconda/bin/aws'
        }
        region = 'us-east-1'
    }

Different queues bound to the same or different Compute environments can be configured according to each process' requirements.


Container Options
=================

As of version ``21.12.1-edge``, the use of the Nextflow :ref:`process-containerOptions` directive is supported to fine control
the properties of the container execution associated with each Batch job.

Not all the standard container options are supported by AWS Batch. These are the options accepted ::

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

Container options must be passed in their long from for "--option value" or short form "-o value", if available.

Few examples::

  containerOptions '--tmpfs /run:rw,noexec,nosuid,size=128 --tmpfs /app:ro,size=64'

  containerOptions '-e MYVAR1 --env MYVAR2=foo2 --env MYVAR3=foo3 --memory-swap 3240000 --memory-swappiness 20 --shm-size 16000000'

  containerOptions '--ulimit nofile=1280:2560 --ulimit nproc=16:32 --privileged'

Check the `AWS doc <https://docs.aws.amazon.com/batch/latest/APIReference/API_ContainerProperties.html>`_ for further details.


Custom AMI
==========

There are several reasons why you might need to create your own `AMI (Amazon Machine Image) <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html>`_
to use in your Compute environments. Typically:

- You do not want to modify your existing Docker images and prefer to install the CLI tool on the hosting environment

- The existing AMI (selected from the marketplace) does not have Docker installed

- You need to attach a larger storage to your EC2 instance (the default ECS instance AMI has only a 30G storage volume which may not be enough for most data analysis pipelines)

- You need to install additional software, not available in the Docker image used to execute the job

Create your custom AMI
----------------------

In the EC2 Dashboard, click the `Launch Instance` button, then choose `AWS Marketplace` in the left pane and enter
`ECS` in the search box. In result list select `Amazon ECS-Optimized Amazon Linux 2 AMI`, then continue as usual to
configure and launch the instance.

.. note:: The selected instance has a bootstrap volume of 8GB and a second EBS volume 30G for computation which is
  hardly enough for real world genomic workloads. Make sure to specify an amount of storage in the second volume
  large enough for the needs of your pipeline execution.

When the instance is running, SSH into it (or connect with the Session Manager service), install the AWS CLI tool
or any other tool that may be required (see next sections).

Once done that, create a new AMI by using the *Create Image* option in the EC2 Dashboard or the AWS command line tool.

The new AMI ID needs to be specified when creating the Batch Compute Environment.

.. warning:: Any installation must be completed on the EC2 instance BEFORE creating the AMI.

.. _aws-cli:

AWS CLI installation
--------------------

.. warning:: The `AWS CLI tool <https://aws.amazon.com/cli>`_ must to be installed in your custom AMI
  by using a self-contained package manager such as `Conda <https://conda.io>`_.

The reason is that when the AWS CLI tool executes using Conda it will use the version of python supplied by Conda.
If you don't use Conda and install the AWS CLI using something like `pip <https://pypi.org/project/pip/>`_ the ``aws``
command will attempt to run using the version of python found in the running container which won't be able to find
the necessary dependencies.

The following snippet shows how to install AWS CLI with `Miniconda <https://conda.io/miniconda.html>`_ in the home folder::

    cd $HOME
    sudo yum install -y bzip2 wget
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
    $HOME/miniconda/bin/conda install -c conda-forge -y awscli
    rm Miniconda3-latest-Linux-x86_64.sh

When complete, verify that the AWS CLI package works correctly::

    $ ./miniconda/bin/aws --version
    aws-cli/1.19.79 Python/3.8.5 Linux/4.14.231-173.361.amzn2.x86_64 botocore/1.20.79

.. note:: The ``aws`` tool will be placed in a directory named ``bin`` in the main installation folder.
  Modifying this directory structure, after the installation, will cause the tool to not work properly.

To configure Nextflow to use this installation, specify the ``cliPath`` parameter in the :ref:`AWS Batch<config-aws-batch>`
configuration as shown below::

    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'

Replace the path above with the one matching the location where ``aws`` tool is installed in your AMI.

.. warning:: The grandparent directory of the ``aws`` tool will be mounted into the container at the same path as the host,
  e.g. ``/home/ec2-user/miniconda``, which will shadow existing files in the container.
  Ensure you use a path that is not already present in the container.

.. note:: Using a version of Nextflow prior 19.07.x the config setting `executor.awscli` should be used
  instead of `aws.batch.cliPath`.

Docker installation
---------------------------------------
Docker is required by Nextflow to execute tasks on AWS Batch. `Amazon ECS-Optimized Amazon Linux 2 AMI` has Docker installed,
however if you create your AMI starting from a different AMI that does not have Docker installed, you need to do it manually.

The following snippet shows how to install Docker on an Amazon EC2 instance::

    sudo yum update -y
    sudo amazon-linux-extras install docker
    sudo yum install docker
    sudo service docker start

Then, add the ``ec2-user`` to the docker group so you can execute Docker commands without using ``sudo``::

    sudo usermod -a -G docker ec2-user

You may have to reboot your instance to provide permissions for the ``ec2-user`` to access the Docker daemon. This has
to be done BEFORE creating the AMI from the current EC2 instance.

Amazon ECS container agent installation
---------------------------------------
The `ECS container agent <https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ECS_agent.html>`_ is a component
of Amazon Elastic Container Service (Amazon ECS) and is responsible for managing containers on behalf of Amazon ECS.
AWS Batch uses Amazon ECS to execute containerized jobs and therefore requires the agent to be installed on compute
resources within your Compute environments.

The ECS container agent is included in the `Amazon ECS-Optimized Amazon Linux 2 AMI`, but if you select a different AMI
you can also install it on any EC2 instance that supports the Amazon ECS specification.

To install the agent, follow these steps::

    sudo amazon-linux-extras disable docker
    sudo amazon-linux-extras install -y ecs
    sudo systemctl enable --now ecs

To test the installation::

    curl -s http://localhost:51678/v1/metadata | python -mjson.tool (test)

.. note:: The ``AmazonEC2ContainerServiceforEC2Role`` policy must be attached to the instance role in order to be able to
    connect the EC2 instance created by the Compute Environment to the ECS container.


Jobs & Execution
================

Custom job definition
---------------------

Nextflow automatically creates the Batch `Job definitions <http://docs.aws.amazon.com/batch/latest/userguide/job_definitions.html>`_
needed to execute your pipeline processes. Therefore it's not required to define them before running your workflow.

However you may still need to specify a custom `Job Definition` to fine control the configuration settings
of a specific job e.g. to define custom mount paths or other Batch Job special settings.

To do that first create a *Job Definition* in the AWS Console (or with other means). Note the name of the *Job Definition*
you created. You can then associate a process execution with this *Job definition* by using the :ref:`process-container`
directive and specifing, in place of the container image name, the Job definition name prefixed by the
``job-definition://`` string, as shown below::

  process.container = 'job-definition://your-job-definition-name'

Pipeline execution
------------------

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads.

Pipeline input data can be stored either locally or in a `S3 <https://aws.amazon.com/s3/>`_ bucket.
The pipeline execution must specifies a AWS Storage bucket where jobs intermediate results are stored with the
``-bucket-dir`` command line options. For example::

  nextflow run my-pipeline -bucket-dir s3://my-bucket/some/path

.. warning::
  The bucket path should include at least a top level directory name e.g. use ``s3://my-bucket/work``
  not just ``s3://my-bucket``.

Hybrid workloads
----------------

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment
of hybrid workloads in which some jobs are executed in the local computer or local computing cluster and
some jobs are offloaded to AWS Batch service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the AWS Batch :ref:`configuration <aws-batch-config>` only to a subset of processes in your workflow.
For example::

  aws {
      region = 'eu-west-1'
      batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
      }
  }

  process {
      withLabel: bigTask {
        executor = 'awsbatch'
        queue = 'my-batch-queue'
        container = 'my/image:tag'
    }
  }

The above configuration snippet will deploy the execution with AWS Batch only for processes annotated
with the :ref:`process-label` ``bigTask``, the remaining process with run in the local computer.

Volume mounts
-------------

User provided container volume mounts can be provided as shown below::

  aws {
    region = 'eu-west-1'
    batch {
        volumes = '/tmp'
    }
  }

Multiple volumes can be specified using a comma separated paths. The usual Docker volume mount syntax
can be used to specify complex volumes for which the container paths is different from the host paths
or to specify *read-only* option. For example::

  aws {
    region = 'eu-west-1'
    batch {
        volumes = ['/tmp', '/host/path:/mnt/path:ro']
    }
  }

The above snippet defines two volume mounts the jobs executed in your pipeline. The first mounting the
host path ``/tmp`` in the same path in the container and using *read-write* access mode. The second
mounts the path ``/host/path`` in the host environment to the ``/mnt/path`` in the container using the
*read-only* access mode.

.. note:: This feature requires Nextflow version 19.07.x or later.

Troubleshooting
---------------

**Problem**: The Pipeline execution terminates with an AWS error message similar to the one shown below::

    JobQueue <your queue> not found

Make sure you have defined a AWS region in the Nextflow configuration file and it matches the region
in which your Batch environment has been created.

**Problem**: A process execution fails reporting the following error message::

  Process <your task> terminated for an unknown reason -- Likely it has been terminated by the external system

This may happen when Batch is unable to execute the process script. A common cause of this problem is that the
Docker container image you have specified uses a non standard `entrypoint <https://docs.docker.com/engine/reference/builder/#entrypoint>`_
which does not allow the execution of the Bash launcher script required by Nextflow to run the job.

This may also happen if the AWS CLI doesn't run correctly.

Other places to check for error information:

- The ``.nextflow.log`` file.
- The Job execution log in the AWS Batch dashboard.
- The `CloudWatch <https://aws.amazon.com/cloudwatch/>`_ logs found in the ``/aws/batch/job`` log group.

**Problem**: A process execution is stalled in the ``RUNNABLE`` status and the pipeline output is similar to the one below::

    executor >  awsbatch (1)
    process > <your process> (1) [  0%] 0 of ....

It may happen that the pipeline execution hangs indefinitely because one of the jobs is held in the queue and never gets
executed. In AWS Console, the queue reports the job as ``RUNNABLE`` but it never moves from there.

There are multiple reasons why this can happen. They are mainly related to the Compute Environment workload/configuration,
the docker service or container configuration, network status, etc.

This `AWS page <https://aws.amazon.com/premiumsupport/knowledge-center/batch-job-stuck-runnable-status/>`_ provides several
resolutions and tips to investigate and work around the issue.


Advanced configuration
======================

Read :ref:`AWS Batch configuration<config-aws-batch>` section to learn more about advanced Batch configuration options.
