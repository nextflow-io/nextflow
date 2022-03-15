.. _executor-page:

*********
Executors
*********

In the Nextflow framework architecture, the `executor` is the component that determines the system where a pipeline
process is run and supervises its execution.

The `executor` provides an abstraction between the pipeline processes and the underlying execution system. This
allows you to write the pipeline functional logic independently from the actual processing platform.

In other words you can write your pipeline script once and have it running on your computer, a cluster resource manager
or the cloud by simply changing the executor definition in the Nextflow configuration file.


.. _awsbatch-executor:

AWS Batch
=========

Nextflow supports `AWS Batch <https://aws.amazon.com/batch/>`_ service which allows submitting jobs in the cloud
without having to spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks,
which makes deploying pipelines much simpler.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'awsbatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads. Moreover a S3 bucket must be used as pipeline work directory.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`
* :ref:`process-queue`

See the :ref:`AWS Batch<aws-batch>` page for further configuration details.


.. _azurebatch-executor:

Azure Batch
===========

Nextflow supports `Azure Batch <https://azure.microsoft.com/en-us/services/batch/>`_ service which allows submitting jobs in the cloud
without having to spin out and manage a cluster of virtual machines. Azure Batch uses Docker containers to run tasks,
which makes deploying pipelines much simpler.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'azurebatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a cloud virtual machine. The latter is suggested for heavy or long
running workloads. Moreover a Azure Blob storage container must be used as pipeline work directory.

See the :ref:`Azure Batch <azure-batch>` page for further configuration details.


.. _ga4ghtes-executor:

GA4GH TES
=========

.. warning:: This is an experimental feature and it may change in a future release. It requires Nextflow
  version 0.31.0 or later.

The `Task Execution Schema <https://github.com/ga4gh/task-execution-schemas>`_ (TES) project
by the `GA4GH <https://www.ga4gh.org>`_ standardisation initiative is an effort to define a
standardized schema and API for describing batch execution tasks in portable manner.

Nextflow includes an experimental support for the TES API providing a ``tes`` executor which allows
the submission of workflow tasks to a remote execution back-end exposing a TES API endpoint.

To use this feature define the following variables in the workflow launching environment::

    export NXF_MODE=ga4gh
    export NXF_EXECUTOR=tes
    export NXF_EXECUTOR_TES_ENDPOINT='http://back.end.com'

It is important that the endpoint is specified without the trailing slash; otherwise, the resulting URLs will be not
normalized and the requests to TES will fail.

Then you will be able to run your workflow over TES using the usual Nextflow command line. Be sure to specify the Docker
image to use, i.e.::

    nextflow run rnaseq-nf -with-docker alpine

.. note:: If the variable ``NXF_EXECUTOR_TES_ENDPOINT`` is omitted the default endpoint is ``http://localhost:8000``.

.. tip:: You can use a local `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ server using the following launch
  command line::

  ./funnel server --Server.HTTPPort 8000 --LocalStorage.AllowedDirs $HOME run

  (tested with version 0.8.0 on macOS)

.. warning:: Make sure the TES back-end can access the workflow work directory when
  data is exchanged using a local or shared file system.

**Known Limitations**

* Automatic deployment of workflow scripts in the `bin` folder is not supported.
* Process output directories are not supported. For details see `#76 <https://github.com/ga4gh/task-execution-schemas/issues/76>`_.
* Glob patterns in process output declarations are not supported. For details see `#77 <https://github.com/ga4gh/task-execution-schemas/issues/77>`_.


.. _google-lifesciences-executor:

Google Life Sciences
====================

`Google Cloud Life Sciences <https://cloud.google.com/life-sciences>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for the Life Sciences API which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as pipelines (it requires Nextflow 20.01.0 or later).

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor set the property ``process.executor = 'google-lifesciences'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-machineType`
* :ref:`process-memory`

See the :ref:`Google Life Sciences <google-lifesciences>` page for further configuration details.


.. _htcondor-executor:

HTCondor
========

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

The ``condor`` executor allows you to run your pipeline script by using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``condor_submit`` command.

Being so, the pipeline must be launched from a node where the ``condor_submit`` command is available, that is, in a
common usage scenario, the cluster `head` node.

.. note::
  The HTCondor executor for Nextflow does not support at this time the HTCondor ability to transfer input/output data to
  the corresponding job computing node. Therefore the data needs to be made accessible to the computing nodes using
  a shared file system directory from where the Nextflow workflow has to be executed (or specified via the ``-w`` option).

To enable the HTCondor executor simply set ``process.executor = 'condor'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`
* :ref:`process-time`


.. _ignite-executor:

Ignite
======

.. danger::
  This feature has been phased out and is no longer supported as of version 22.01.x.

The ``ignite`` executor allows you to run a pipeline on an `Apache Ignite <https://ignite.apache.org/>`_ cluster.

To enable this executor set ``process.executor = 'ignite'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`

See the :ref:`ignite-page` page to learn how to configure Nextflow to deploy and run an
Ignite cluster in your infrastructure.


.. _k8s-executor:

Kubernetes
==========

The ``k8s`` executor allows you to run a pipeline on a `Kubernetes <http://kubernetes.io/>`_ cluster.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-pod`

See the :ref:`Kubernetes <k8s-page>` page to learn how to set up a Kubernetes cluster for running Nextflow pipelines.


.. _local-executor:

Local
=====

The ``local`` executor is used by default. It runs the pipeline processes in the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads` and by taking advantage of multi-cores
architecture provided by the CPU.

In a common usage scenario, the `local` executor can be useful to develop and test your pipeline script in your computer,
switching to a cluster facility when you need to run it on production data.


.. _lsf-executor:

LSF
===

The ``lsf`` executor allows you to run your pipeline script by using a `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``bsub`` command.

Being so, the pipeline must be launched from a node where the ``bsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the LSF executor simply set ``process.executor = 'lsf'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. note::

    LSF supports both *per-core* and *per-job* memory limit. Nextflow assumes that LSF works in the
    *per-core* memory limits mode, thus it divides the requested :ref:`process-memory` by the number of requested :ref:`process-cpus`.

    This is not required when LSF is configured to work in *per-job* memory limit mode. You will need to specified that
    adding the option ``perJobMemLimit`` in :ref:`config-executor` in the Nextflow configuration file.

    See also the `Platform LSF documentation <https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita>`_.


.. _moab-executor:

Moab
====

The ``moab`` executor allows you to run your pipeline script by using the
`Moab <https://en.wikipedia.org/wiki/Moab_Cluster_Suite>`_ resource manager by
`Adaptive Computing <http://www.adaptivecomputing.com/>`_.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``msub`` command provided
by the resource manager.

Being so, the pipeline must be launched from a node where the ``msub`` command is available, that is, in a common usage
scenario, the compute cluster `login` node.

To enable the `Moab` executor simply set ``process.executor = 'moab'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _nqsii-executor:

NQSII
=====

The ``nsqii`` executor allows you to run your pipeline script by using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the NQSII executor simply set ``process.executor = 'nqsii'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _oar-executor:

OAR
===

The ``oar`` executor allows you to run your pipeline script by using the `OAR <https://oar.imag.fr>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``oarsub`` command.

Being so, the pipeline must be launched from a node where the ``oarsub`` command is available, that is, in a common usage scenario, the cluster `head` node.

To enable the OAR executor simply set ``process.executor = 'oar'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

**Known Limitations**

* ``clusterOptions`` should be given, if more than one, semicolon separated. It ensures the `OAR` batch script to be accurately formatted::

    clusterOptions = '-t besteffort;--project myproject'


.. _pbs-executor:

PBS/Torque
==========

The ``pbs`` executor allows you to run your pipeline script by using a resource manager belonging to the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ family of batch schedulers.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the PBS executor simply set ``process.executor = 'pbs'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _pbspro-executor:

PBS Pro
=======

The ``pbspro`` executor allows you to run your pipeline script by using the `PBS Pro <https://www.pbspro.org/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``qsub`` command provided
by the scheduler.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `login` node.

To enable the PBS Pro executor simply set ``process.executor = 'pbspro'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _sge-executor:

SGE
===

The ``sge`` executor allows you to run your pipeline script by using a `Sun Grid Engine <http://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_
cluster or a compatible platform (`Open Grid Engine <http://gridscheduler.sourceforge.net/>`_, `Univa Grid Engine <http://www.univa.com/products/grid-engine.php>`_, etc).

Nextflow manages each process as a separate grid job that is submitted to the cluster by using the ``qsub`` command.

Being so, the pipeline must be launched from a node where the ``qsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SGE executor simply set ``process.executor = 'sge'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-penv`
* :ref:`process-queue`
* :ref:`process-time`


.. _slurm-executor:

SLURM
=====

The ``slurm`` executor allows you to run your pipeline script by using the `SLURM <https://slurm.schedmd.com/documentation.html>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``sbatch`` command.

Being so, the pipeline must be launched from a node where the ``sbatch`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SLURM executor simply set ``process.executor = 'slurm'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. note:: SLURM `partitions` can be considered jobs queues. Nextflow allows you to set partitions by using the above ``queue``
    directive.

.. tip:: Nextflow does not provide a direct support for SLURM multi-clusters feature. If you need to
  submit workflow executions to a cluster that is not the current one, specify it setting the
  ``SLURM_CLUSTERS`` variable in the launching environment.
