.. _awscloud-page:

************
Amazon Cloud
************

AWS credentials
---------------

Nextflow will use the AWS credentials defined in your environment, using the standard AWS variables shown below:

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


.. _awscloud-batch:

AWS Batch
=========

.. note::
    Requires Nextflow version `0.26.0` or later.

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
    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'

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
from the host image specify the ``cliPath`` parameter in the :ref:`AWS Batch<config-aws-batch>`
configuration as shown below::

    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'

Replace the path above with the one matching the location where ``aws`` tool is installed in your AMI.

.. note:: Using a version of Nextflow prior 19.07.x the config setting `executor.awscli` should be used
  instead of `aws.batch.cliPath`.

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
of hybrid workloads in which some jobs are execute in the local computer or local computing cluster and
some jobs are offloaded to AWS Batch service.

To enable this feature use one or more :ref:`config-process-selectors` in your Nextflow configuration file to apply
the AWS Batch :ref:`configuration <awscloud-batch-config>` only to a subset of processes in your workflow.
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

Advanced configuration
----------------------

Read :ref:`AWS Batch configuration<config-aws-batch>` section to learn more about advanced Batch configuration options.
