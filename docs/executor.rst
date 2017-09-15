.. _executor-page:

***********
Executors
***********

In the Nextflow framework architecture, the `executor` is the component that determines the system where a pipeline
process is run and supervises its execution.

The `executor` provides an abstraction between the pipeline processes and the underlying execution system. This
allows you to write the pipeline functional logic independently from the actual processing platform.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager
or the cloud by simply changing the executor definition in the Nextflow configuration file.

.. _local-executor:

Local
=====

The `local` executor is used by default. It runs the pipeline processes in the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads` and by taking advantage of multi-cores
architecture provided by the CPU.

In a common usage scenario, the `local` executor can be useful to develop and test your pipeline script in your computer,
switching to a cluster facility when you need to run it on production data.


.. _sge-executor:

SGE
===

The `SGE` executor allows you to run your pipeline script by using a `Sun Grid Engine <http://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_
cluster or a compatible platform (`Open Grid Engine <http://gridscheduler.sourceforge.net/>`_, `Univa Grid Engine <http://www.univa.com/products/grid-engine.php>`_, etc).

Nextflow manages each process as a separate grid job that is submitted to the cluster by using the ``qsub`` command.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SGE executor simply set to ``process.executor`` property to ``sge`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-memory`
* :ref:`process-penv`
* :ref:`process-time`
* :ref:`process-clusterOptions`

.. _lsf-executor:

LSF
===

The `LSF` executor allows you to run your pipeline script by using a `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``bsub`` command.

Being so, the pipeline must be launched from a node where the ``bsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the LSF executor simply set to ``process.executor`` property to ``lsf`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. note::

    LSF supports both *per-core* and *per-job* memory limit. Nextflow assumes that LSF works in the
    *per-core* memory limits mode, thus it divides the requested :ref:`process-memory` by the number of requested :ref:`process-cpus`.

    This is not required when LSF is configured to work in *per-job* memory limit mode. You will need to specified that
    adding the option ``perJobMemLimit`` in :ref:`config-executor` in the Nextflow configuration file.

    See also the `Platform LSF documentation <https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita>`_.


.. _slurm-executor:

SLURM
=====


The `SLURM` executor allows you to run your pipeline script by using the `SLURM <https://computing.llnl.gov/linux/slurm/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``sbatch`` command.

Being so, the pipeline must be launched from a node where the ``sbatch`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SLURM executor simply set to ``process.executor`` property to ``slurm`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. note:: SLURM `partitions` can be considered jobs queues. Nextflow allows to set partitions by using the above ``queue``
    directive.

.. _pbs-executor:

PBS/Torque
==========

The `PBS` executor allows you to run your pipeline script by using a resource manager belonging to the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ family of batch schedulers.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the PBS executor simply set the property ``process.executor = 'pbs'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _nqsii-executor:

NQSII
=====

The `NQSII` executor allows you to run your pipeline script by using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the NQSII executor simply set the property ``process.executor = 'nqsii'`` in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`

.. _condor-executor:

HTCondor
========

The `HTCondor` executor allows you to run your pipeline script by using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ resource manager.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``condor_submit`` command.

Being so, the pipeline must be launched from a node where the ``condor_submit`` command is available, that is, in a
common usage scenario, the cluster `head` node.

To enable the HTCondor executor simply set to ``process.executor`` property to ``condor`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-disk`
* :ref:`process-clusterOptions`


.. _drmaa-executor:

DRMAA
=====

The `DRMAA` executor allows you to execute a Nextflow pipeline by using a grid engine that implements the
Java binding for the `DRMAA <http://www.drmaa.org>`_ interface api (version 1).

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

In order to be able to use this executor you will need to access the DRMAA libraries provided by your cluster vendor.
Commonly these files are named ``drmaa.jar`` and ``libdrmaa.so`` and they are located in the cluster installation lib folder.
Ask to your IT administrator how to find these files.

To enable the PBS executor you will need to set the property ``process.executor='drmaa'`` in the ``nextflow.config`` file,
moreover you will need to specify the ``drmaa.jar`` library path on the Nextflow command line by using the ``-with-drmaa``
option. For example::

  nextflow run <your pipeline> -with-drmaa /some/path/drmaa.jar


Alternatively, instead of specifying the DRMAA library on the command line, you may want to use the environment variable
``NXF_DRMAA`` to define it.

.. tip:: If you get the following error message::

      ERROR: java.lang.UnsatisfiedLinkError: no drmaa in java.library.path

    It means that Nextflow is unable to find ``libdrmaa.so`` file. The most common solution is
    to include the path where this file is located in the ``LD_LIBRARY_PATH`` environment variable.


The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`


.. _dnanexus-executor:

DNAnexus
========

The `DNAnexus` executor allows you to run your pipeline in the `DNAnexus <http://dnanexus.com/>`_ cloud platform.

Nextflow pipeline to be executed in the DNAnexus platform need to be packaged as DNAnexus app. Read how bundle and
deploy Nextflow apps in the :ref:`dnanexus-page` section.

The `dnanexus` executor allows your script to submit pipeline's processes in the DNAnexus cloud platform as separate jobs.
It has to be specified in the configuration file or on the program command line options, as shown below::

    process.executor = 'dnanexus'


Property 'instanceType'
------------------------

The ``instanceType`` configuration property allows you to specify the instance type to be used by a process when
executing the required job. For example::

     process.executor = 'dnanexus'
     process.instanceType = 'dx_m1.xlarge'


The list of instance types, that can be used for this property, is available in the `Run Specification
<https://wiki.dnanexus.com/API-Specification-v1.0.0/IO-and-Run-Specifications#Run-Specification>`_ page.

.. _ignite-executor:

Ignite
======

The `Ignite` executor allows you to run a pipeline by using the `Apache Ignite <https://ignite.apache.org/>`_ clustering
technology that is embedded with the Nextflow runtime.

To enable this executor set the property ``process.executor = 'ignite'`` in the ``nextflow.config`` file.

The amount of resources requested by each task submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`

Read the :ref:`ignite-page` section in this documentation to learn how to configure Nextflow to deploy and run an
Ignite cluster in your infrastructure.


Kubernetes
==========

Nextflow provides an experimental support for `Kubernetes <http://kubernetes.io/>`_ clustering technology. It allows
you to deploy and transparently run a Nextflow pipeline in a Kubernetes cluster.

Nextflow manages each process as a separate `pod` that is submitted to the cluster by using the ``kubectl`` command.
Being so, the pipeline must be launched from a node where the ``kubectl`` command is available.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover the pipeline must be executed in a shared file system folder
accessible by all Kubernetes cluster nodes.

To enable this executor set the property ``process.executor = 'k8s'`` in the ``nextflow.config`` file.

The following directives can be used to define the amount of computing resources needed and the container(s) to use:

* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-container`


AWS Batch (Beta)
================

Nextflow supports `AWS Batch <https://aws.amazon.com/batch/>`_ service which allows submitting jobs in the cloud without having to
spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks, which makes deploying pipelines much simpler.


Installation 
------------

In order to use AWS Batch, it is necessary to clone the Nextflow repository, move to the ``aws-batch`` branch and compile Nextflow::

    git clone https://github.com/nextflow-io/nextflow.git
    git checkout aws-batch
    make compile


Configuration
-------------

To submit jobs on AWS Batch it is necessary to specify ``aws-batch`` as executor in the ``nextflow.config`` file.

Once done that a ``queue`` parameter needs to specified, either in ``nextflow.config`` or within each process of the pipeline::

    process {
        queue = 'my-queue/my-job-definition'
    }

Finally to run Nextflow with AWS Batch, an S3 bucket needs to be used as a working directory::

    nextflow run pipeline.nf -w s3://my-bucket/my-prefix

The official AWS Batch documentation provides all the instructions to create `computing environments <http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html>`_, `queues <http://docs.aws.amazon.com/batch/latest/userguide/job_queues.html>`_ and `job definitions <http://docs.aws.amazon.com/batch/latest/userguide/job_definitions.html>`_.

AWS Batch uses Amazon ECS service to spin up clusters of Docker containers, so each job definiton needs to be linked to an existing container that will be used to run
the tasks. Docker container images need to be available in a repository that can be reached by the ECS instances (so either on Docker Hub or on AWS ECR for example).

Nextflow needs also to use the ``aws cli`` tool in order to copy the runnable scripts inside the container and perform the staging and unstaging of files with S3.
The aws cli can either be already present or added to existing containers or can be installed directly on a custom AMI that will be used by Batch (see next section).


Custom AMI
----------

AWS Batch uses the default ECS instance AMI, which has only a 22 GB volume for Docker container, so a custom AMI needs to be prepared or used in order to fully utilize the service to run bioinformatics pipelines with Nextflow.

In order to do that, follow the instruction from AWS documentation to create a `custom AMI <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/creating-an-ami-ebs.html>`_. It is highly reccommended to create also an EBS volume to be used as temporary scratch space for Nextflow tasks.

So let's assume that the custom AMI has a EBS volume attached and mounted as a partition called ``/scratch``, then in the ``job definition`` for AWS Batch this volume on the host instance needs to be mounted as ``/tmp`` on the container. In this way the Nextflow task that will run inside the container will directly use the scratch space on the host for all the intermediate files. 

If the ``aws cli`` tool needs to be installed outside the Docker containers, during the custom AMI creation it is possible to put a self-contained python environment (like virtualenv or MiniConda) containing the ``aws cli``. For example to do that with Miniconda::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /scratch/miniconda
    /scratch/miniconda/bin/conda install -c conda-forge awscli

Once done that, in order to work from the container, the ``aws cli`` executable needs to be modified to reflect the path of the mount point in the container. So, as suggested, if the host ``/scratch`` volume will be mounted as ``/tmp`` in the container, then the shebang of the executable ``/scratch/miniconda/bin/aws`` needs to be changed like this::

   #!/tmp/miniconda/bin/python3.6

By default Nextflow will assume the ``aws cli`` is directly available in the container. To use an installation from the host image it is possible to specify a ``aws_cli`` parameter in the Nextflow process configuration::

    process {
        aws_cli = '/tmp/miniconda/bin/aws'
    }


