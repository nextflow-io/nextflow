.. _awscloud-page:

************
Amazon Cloud
************

Nextflow provides out of the box support for the Amazon AWS cloud allowing you to setup a computing cluster,
deploy it and run your pipeline in the AWS infrastructure in a few commands.


Configuration
=============

Cloud configuration attributes are provided in the ``nextflow.config`` file as shown in the example below::

    cloud {
        imageId = 'ami-4b7daa32'
        instanceType = 'm4.xlarge'
        subnetId = 'subnet-05222a43'
    }

The above attributes define the virtual machine ID and type to be used and the VPC subnet ID to be applied
in you cluster. Replace these values with the ones of your choice.

Nextflow only requires a Linux image that provides support for `Cloud-init <http://cloudinit.readthedocs.io/>`_
bootstrapping mechanism and includes a Java runtime (version 8) and a Docker engine (version 1.11 or higher).

For your convenience the following pre-configured Amazon Linux AMI is available in the *EU Ireland* region:
``ami-4b7daa32``.


AWS credentials
---------------

Nextflow will use the AWS credentials defined in your environment, using the standard AWS variables shown below:

    * ``AWS_ACCESS_KEY_ID``
    * ``AWS_SECRET_ACCESS_KEY``
    * ``AWS_DEFAULT_REGION``


Alternatively AWS credentials can be specified in the Nextflow configuration file.
See :ref:`AWS configuration<config-aws>` for more details.

.. note:: Credentials can also be provided by using an IAM Instance Role. The benefit of this approach is that
  it spares you from managing/distributing AWS keys explicitly.
  Read the `IAM Roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_ documentation
  and `this blog post <https://aws.amazon.com/blogs/security/granting-permission-to-launch-ec2-instances-with-iam-roles-passrole-permission/>`_ for more details.

User & SSH key
--------------

By default Nextflow creates in each EC2 instance a user with the same name as the one in your local computer and install
the SSH public key available at the path ``$HOME/.ssh/id_rsa.pub``. A different user/key can be specified as shown below::

    cloud {
        userName = 'the-user-name'
        keyFile = '/path/to/ssh/key.pub'
    }

If you want to use a *key-pair* defined in your AWS account and the default user configured in the AMI, specify the
attribute ``keyName`` in place of ``keyFile`` and the name of the existing user specifying the ``userName`` attribute.


Storage
-------

The following storage types can be defined in the cloud instances: *boot*, *instance* and *shared*.

Boot storage
------------

You can set the size of the `root device volume <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/RootDeviceStorage.html>`_
by specifying the attribute ``bootStorageSize``. For example::

    cloud {
        imageId = 'ami-xxx'
        bootStorageSize = '10 GB'
    }


Instance storage
----------------

Amazon instances can provide one or more `ephemeral storage volumes <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html>`_,
depending the instance type chosen. This storage can be made available by using the ``instanceStorageMount``
configuration attributes, as shown below::

    cloud {
        imageId = 'ami-xxx'
        instanceStorageMount = '/mnt/scratch'
    }


The mount path can be any of your choice.

.. note:: When the selected instance provides more than one ephemeral storage volume, Nextflow automatically groups all
  of them together in a single logical volume and mounts it to the specified path. Therefore the resulting instance
  storage size is equals to the sum of the sizes of all ephemeral volumes provided by the actual instance
  (this feature requires Nextflow version 0.27.0 or higher).

If you want to mount a specific instance storage volume, specify the corresponding device name by using
the ``instanceStorageDevice`` setting in the above configuration. See the Amazon documentation for details on
`EC2 Instance storage <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html>`_ and
`devices naming <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/device_naming.html>`_.


Shared file system
------------------

A shared file system can easily be made available in your cloud cluster by using the `Amazon EFS <https://aws.amazon.com/efs/>`_
storage. You will only need to specify in the cloud configuration the EFS file system ID and optionally the
mount path. For example::

    cloud {
        imageId = 'ami-xxx'
        sharedStorageId = 'fs-1803efd1'
        sharedStorageMount = '/mnt/shared'
    }

.. note:: When the attribute ``sharedStorageMount`` is omitted the path ``/mnf/efs`` is used by default.


Cluster deployment
==================

Once defined the configuration settings in the ``nextflow.config`` file you can create the cloud cluster
by using the following command::

    nextflow cloud create my-cluster -c <num-of-nodes>

The string ``my-cluster`` identifies the cluster instance. Replace it with a name of your choice.

Finally replace ``num-of-nodes`` with the actual number of instances that will made-up the cluster.
One node is created as *master*, the remaining as *workers*. If the option ``-c`` is omitted only the *master* node
is created.

.. warning:: You will be charged accordingly the type and the number of instances chosen.


Pipeline execution
==================

Once the cluster initialization is complete, connect to the *master* node using the SSH command which will be displayed by
Nextflow.

.. note:: On MacOS, use the following command to avoid being asked for a pass-phrase even
  you haven't defined one::

    ssh-add -K [private key file]

You can run your Nextflow pipeline as usual, the environment is automatically configured to use the :ref:`Ignite<ignite-page>`
executor. If the Amazon EFS storage is specified in the cloud configuration the Nextflow work directory will
automatically be set in a shared folder in that file system.

The suggested approach is to run your pipeline downloading it from a public repository such
GitHub and to pack the binaries dependencies in a Docker container as described in the
:ref:`Pipeline sharing <sharing-page>` section.

Cluster shutdown
================

When completed shutdown the cluster instances by using the following command::

    nextflow cloud shutdown my-cluster


Cluster auto-scaling
====================

Nextflow integration for AWS cloud provides a native support auto-scaling that allows the computing cluster
to scale-out or scale-down i.e. add or remove computing nodes dynamically at runtime.

This is a critical feature, especially for pipelines crunching not homogeneous dataset, because it allows the
cluster to adapt dynamically to the actual workload computing resources need as they change over the time.

Cluster auto-scaling is enabled by adding the ``autoscale`` option group in the configuration file as shown below::

    cloud {
        imageId = 'ami-xxx'
        autoscale {
            enabled = true
            maxInstances = 10
        }
    }


The above example enables automatic cluster scale-out i.e. new instances are automatically launched and added to the
cluster when tasks remain too long in wait status because there aren't enough computing resources available. The
``maxInstances`` attribute defines the upper limit to which the cluster can grow.

By default unused instances are not removed when are not utilised. If you want to enable automatic cluster scale-down
specify the ``terminateWhenIdle`` attribute in the ``autoscale`` configuration group.

It is also possible to define a different AMI image ID, type and spot price for instances launched by the Nextflow autoscaler.
For example::

    cloud {
        imageId = 'ami-xxx'
        instanceType = 'm4.large'

        autoscale {
            enabled = true
            spotPrice = 0.15
            minInstances = 5
            maxInstances = 10
            imageId = 'ami-yyy'
            instanceType = 'm4.4xlarge'
            terminateWhenIdle = true
        }
    }

By doing that it's is possible to create a cluster with a single node i.e. the master node. Then the autoscaler will
automatically add the missing instances, up to the number defined by the ``minInstances`` attributes. These will have a
different image and type from the master node and will be launched a *spot instances* because the ``spotPrice``
attribute has been specified.


Spot prices
===========

Nextflow includes an handy command to list the current price of EC2 spot instances. Simply type the following
command in your shell terminal::

    nextflow cloud spot-prices

It will print the current spot price for all available instances type, similar to the example below::

    TYPE        PRICE  PRICE/CPU ZONE       DESCRIPTION             CPUS   MEMORY DISK
    t1.micro    0.0044    0.0044 eu-west-1c Linux/UNIX                 1 627.7 MB -
    m4.4xlarge  0.1153    0.0072 eu-west-1a Linux/UNIX (Amazon VPC)   16    64 GB -
    m4.10xlarge 0.2952    0.0074 eu-west-1b Linux/UNIX (Amazon VPC)   40   160 GB -
    m4.large    0.0155    0.0077 eu-west-1b Linux/UNIX (Amazon VPC)    2     8 GB -
    m4.2xlarge  0.0612    0.0077 eu-west-1a Linux/UNIX (Amazon VPC)    8    32 GB -
    m4.xlarge   0.0312    0.0078 eu-west-1a Linux/UNIX (Amazon VPC)    4    16 GB -
    c4.8xlarge  0.3406    0.0095 eu-west-1c Linux/UNIX (Amazon VPC)   36    60 GB -
    m1.xlarge   0.0402    0.0100 eu-west-1b Linux/UNIX                 4    15 GB 4 x 420 GB
    c4.4xlarge  0.1652    0.0103 eu-west-1b Linux/UNIX (Amazon VPC)   16    30 GB -
    c1.xlarge   0.0825    0.0103 eu-west-1a Linux/UNIX                 8     7 GB 4 x 420 GB
    m1.medium   0.0104    0.0104 eu-west-1b Linux/UNIX (Amazon VPC)    1   3.8 GB 1 x 410 GB
    c3.8xlarge  0.3370    0.0105 eu-west-1a Linux/UNIX                32    60 GB 2 x 320 GB
    c3.2xlarge  0.0860    0.0108 eu-west-1c Linux/UNIX                 8    15 GB 2 x 80 GB
    c3.4xlarge  0.1751    0.0109 eu-west-1c Linux/UNIX (Amazon VPC)   16    30 GB 2 x 160 GB
    m3.2xlarge  0.0869    0.0109 eu-west-1c Linux/UNIX (Amazon VPC)    8    30 GB 2 x 80 GB
    r3.large    0.0218    0.0109 eu-west-1c Linux/UNIX                 2  15.2 GB 1 x 32 GB
    :


It's even possible to refine the showed data by specifying a filtering and ordering criteria. For example::

    nextflow cloud spot-prices -sort pricecpu -filter "cpus==4"


It will only print instance types having 4 cpus and sorting them by the best price per cpu.


Advanced configuration
======================

Read :ref:`Cloud configuration<config-cloud>` section to learn more about advanced cloud configuration options.


.. _awscloud-batch:

AWS Batch
=========

.. warning:: This is an experimental feature and it may change in a future release. It requires Nextflow
  version `0.26.0` or later.

`AWS Batch <https://aws.amazon.com/batch/>`_ is a managed computing service that allows the execution of containerised
workloads in the Amazon cloud infrastructure.

Nextflow provides a built-in support for AWS Batch which allows the seamless deployment of a Nextflow pipeline
in the cloud offloading the process executions as Batch jobs.

.. _awscloud-batch-config:

Configuration
-------------

1 - Make sure your pipeline processes specifies one or more Docker containers by using the :ref:`process-container` directive.

2 - Container images need to be published in a Docker registry such as `Docker Hub <https://hub.docker.com/>`_,
`Quay <https://quay.io/>`_ or `ECS Container Registry <https://aws.amazon.com/ecr/>`_ that can be reached
by ECS Batch.

3 - Specify the AWS Batch :ref:`executor<awsbatch-executor>` in the pipeline configuration.

4 - Specify one or more AWS Batch queues for the execution of your pipeline by using the :ref:`process-queue` directive.
Batch queues allow you to bind the execution of a process to a specific computing environment ie. number of CPUs,
type of instances (On-demand or Spot), scaling ability, etc. See the `AWS Batch documentation <http://docs.aws.amazon.com/batch/latest/userguide/create-job-queue.html>`_ to learn
how to setup Batch queues.

5 - Make sure the container image includes the `AWS CLI tool <https://aws.amazon.com/cli>`_ i.e. ``aws``.
Alternatively, it can also be installed in a custom AMI image. See the note below for details.

An example ``nextflow.config`` file is shown below::

    process.executor = 'awsbatch'
    process.queue = 'my-batch-queue'
    process.container = 'quay.io/biocontainers/salmon'
    aws.region = 'eu-west-1'
    
    // NOTE: this setting is only required if the AWS CLI tool is installed in a custom AMI
    executor.awscli = '/home/ec2-user/miniconda/bin/aws'

.. note:: Nextflow requires to access the AWS command line tool (``aws``) from the container in which the job runs
  in order to stage the required input files and to copy back the resulting output files in the
  `S3 storage <https://aws.amazon.com/s3/>`_.

The ``aws`` tool can either be included in container image(s) used by your pipeline execution or
installed in a custom AMI that needs to used in place of the default AMI when configuring the Batch
`Computing environment <http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html>`_.

The latter approach is preferred  because it allows the use of existing Docker images without the need to add
the AWS CLI tool to them. See the sections below to learn how create a custom AMI and install the AWS CLI tool
to it.

.. warning:: AWS Batch uses the default ECS instance AMI, which has only a 22 GB storage volume which may not
  be enough for real world data analysis pipelines.

See the section below to learn how to create a custom AWS Batch custom AMI with a larger storage.

Custom AMI
----------

In the EC2 Dashboard, click the `Launch Instance` button, then choose `AWS Marketplace` in the left pane and enter
`ECS` in the search box. In result list select `Amazon ECS-Optimized Amazon Linux AMI`, then continue as usual to
configure and launch the instance.

.. note:: The selected instance has a bootstrap volume of 8GB and a second EBS volume 22G for computation which is
  hardly enough for real world genomic workloads. Make sure to specify an amount of storage in the second volume
  large enough for the needs of your pipeline execution.

When the instance is running, SSH into it, install the AWS CLI tools as explained below or any other required tool
that may be required.

Also make sure the Docker configuration reflects the amount of storage you have specified when launching the instance
as shown below::

    $ docker info | grep -i data
     Data file:
     Metadata file:
     Data Space Used: 500.2 MB
     Data Space Total: 1.061 TB
     Data Space Available: 1.06 TB
     Metadata Space Used: 733.2 kB
     Metadata Space Total: 1.074 GB
     Metadata Space Available: 1.073 GB

The above example shows the Docker data configuration for a 1000GB EBS data volume. See the `ECS Storage documentation <http://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-ami-storage-config.html>`_
for more details.

.. warning:: The maximum storage size of a single Docker container is by default 10GB, independently the amount of data space available
  in the underlying volume (see `Base device size <https://docs.docker.com/engine/reference/commandline/dockerd/#dmbasesize>`_ for more details).

You can verify the current setting by using this command::

     $ docker info | grep -i base
       Base Device Size: 10.74 GB

If your pipeline needs more storage for a single task execution, you will need to specify the ``dm.basesize`` setting
with a proper value in the ``/etc/sysconfig/docker-storage`` configuration file.
See `here <https://forums.aws.amazon.com/message.jspa?messageID=811761#811761>`_
and `here <https://www.projectatomic.io/blog/2016/03/daemon_option_basedevicesize/>`_ for details.


Once done that, create a new AMI by using the *Create Image* option in the EC2 Dashboard or the AWS command line tool.

The new AMI ID needs to be specified when creating the Batch
`Computing environment <http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html>`_.

.. _aws-cli:

AWS CLI installation
--------------------

.. warning:: The `AWS CLI tool <https://aws.amazon.com/cli>`_ must to be installed in your custom AMI
  by using a self-contained package manager such as `Conda <https://conda.io>`_.

The reason is that when the AWS CLI tool executes using Conda it will use the version of python supplied by Conda.
If you don't use Conda and install the AWS CLI using something like `pip <https://pypi.org/project/pip/>`_ the ``aws``
command will attempt to run using the version of python found in the running container which won't be able to find
the necessary dependencies.

The following snippet shows how to install AWS CLI with `Miniconda <https://conda.io/miniconda.html>`_::

    sudo yum install -y bzip2 wget
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
    $HOME/miniconda/bin/conda install -c conda-forge -y awscli
    rm Miniconda3-latest-Linux-x86_64.sh

When complete verifies that the AWS CLI package works correctly::

    $ ./miniconda/bin/aws --version
    aws-cli/1.11.120 Python/3.6.3 Linux/4.9.43-17.39.amzn1.x86_64 botocore/1.5.83


.. note:: The ``aws`` tool will be placed in a directory named ``bin`` in the main installation folder.
  Modifying this directory structure, after the installation, will cause the tool to not work properly.


By default Nextflow will assume the AWS CLI tool is directly available in the container. To use an installation
from the host image specify the ``awscli`` parameter in the Nextflow :ref:`executor <awsbatch-executor>`
configuration as shown below::

    executor.awscli = '/home/ec2-user/miniconda/bin/aws'

Replace the path above with the one matching the location where ``aws`` tool is installed in your AMI.

Custom job definition
---------------------

Nextflow automatically creates the Batch `Job definitions <http://docs.aws.amazon.com/batch/latest/userguide/job_definitions.html>`_
needed to execute your pipeline processes. Therefore it's not required to define them before run your workflow.

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

Pipeline input data should to be stored in the Input data `S3 storage <https://aws.amazon.com/s3/>`_. In the same
manner the pipeline execution must specifies a S3 bucket where jobs intermediate results are stored with the
``-bucket-dir`` command line options. For example::

  nextflow run my-pipeline -bucket-dir s3://my-bucket/some/path


.. warning::
  The bucket path should include at least a top level directory name e.g. use ``s3://my-bucket/work``
  not just ``s3://my-bucket``. 

Hybrid workloads
----------------

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment
of hybrid workloads in which some jobs are execute in the local computer or local computing cluster and
some jobs are offloaded to AWS Batch service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the AWS Batch :ref:`configuration <awscloud-batch-config>` only to a subset of processes in your workflow.
For example::


  aws {
      region = 'eu-west-1'
  }
  executor {
    awscli = '/home/ec2-user/miniconda/bin/aws'
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




