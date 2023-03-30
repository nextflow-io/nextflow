.. _executor-page:

*********
Executors
*********

In the Nextflow framework architecture, the `executor` is the component that determines the system where a pipeline
process is run and supervises its execution.

The `executor` provides an abstraction between the pipeline processes and the underlying execution system. This
allows you to write the pipeline functional logic independently from the actual processing platform.

In other words, you can write your pipeline script once and have it running on your computer, a cluster resource manager,
or the cloud — simply change the executor definition in the Nextflow configuration file.


.. _awsbatch-executor:

AWS Batch
=========

Nextflow supports the `AWS Batch <https://aws.amazon.com/batch/>`_ service that allows job submission in the cloud
without having to spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks,
which greatly simplifies pipeline deployment.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor, set the property ``process.executor = 'awsbatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer, or an EC2 instance. EC2 is suggested for heavy or long-running workloads. Moreover, an S3 bucket must be used as the pipeline work directory.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

See the :ref:`AWS Batch<aws-batch>` page for further configuration details.


.. _azurebatch-executor:

Azure Batch
===========

Nextflow supports the `Azure Batch <https://azure.microsoft.com/en-us/services/batch/>`_ service that allows job submission in the cloud
without having to spin out and manage a cluster of virtual machines. Azure Batch uses Docker containers to run tasks,
which greatly simplifies pipeline deployment.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor, set the property ``process.executor = 'azurebatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer, or a cloud virtual machine. The cloud VM is suggested for heavy or long-running workloads. Moreover, an Azure Blob storage container must be used as the pipeline work directory.

See the :ref:`Azure Batch <azure-batch>` page for further configuration details.


.. _bridge-executor:

Bridge
======

`Bridge <https://github.com/cea-hpc/bridge>`_ is an abstraction layer to ease batch system and resource manager usage in
heterogeneous HPC environments.

It is open source software that can be installed on top of existing classical job schedulers such as Slurm, LSF, or other
schedulers. Bridge allows you to submit jobs, get information on running jobs, stop jobs, get information on the cluster system, etc.

For more details on how to install the Bridge system, see the `documentation <https://github.com/cea-hpc/bridge>`_.

To enable the Bridge executor, simply set ``process.executor = 'bridge'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _flux-executor:

Flux Framework Executor
=======================

The ``flux`` executor allows you to run your pipeline script using the `Flux Framework <https://flux-framework.org>`_.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``flux mini submit`` command.

To enable the Flux executor, simply set ``process.executor = 'flux'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-queue`
* :ref:`process-time`

Additionally, to have Flux print all output to stderr and stdout, set `flux.terminalOutput` to true.

.. note:: Flux does not support specifying memory. 


.. _ga4ghtes-executor:

GA4GH TES
=========

.. warning:: This is an experimental feature and it may change in future releases. It requires Nextflow
  version 0.31.0 or later.

The `Task Execution Schema <https://github.com/ga4gh/task-execution-schemas>`_ (TES) project
by the `GA4GH <https://www.ga4gh.org>`_ standardization initiative is an effort to define a
standardized schema and API for describing batch execution tasks in a portable manner.

Nextflow includes experimental support for the TES API by providing a ``tes`` executor, which allows
the submission of workflow tasks to a remote execution back-end exposing a TES API endpoint.

To use this feature, define the following variables in the workflow launching environment::

    export NXF_MODE=ga4gh
    export NXF_EXECUTOR=tes
    export NXF_EXECUTOR_TES_ENDPOINT='http://back.end.com'

It is important that the endpoint is specified without the trailing slash; otherwise, the resulting URLs will not be
normalized and the requests to TES will fail.

You will then be able to run your workflow over TES using the usual Nextflow command line. Be sure to specify the Docker
image to use, i.e.::

    nextflow run rnaseq-nf -with-docker alpine

.. note:: If the variable ``NXF_EXECUTOR_TES_ENDPOINT`` is omitted, the default endpoint is ``http://localhost:8000``.

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

.. _google-batch-executor:

Google Cloud Batch
===================

`Google Cloud Batch <https://cloud.google.com/batch>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for the Batch API that allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as pipelines (requires Nextflow ``22.07.1-edge`` or later).

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover, the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor, set the property ``process.executor = 'google-batch'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-container`
* :ref:`process-containerOptions`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-machineType`
* :ref:`process-memory`
* :ref:`process-time`
* :ref:`process-resourcelabels`

See the :ref:`Google Cloud Batch <google-batch>` page for further configuration details.

.. _google-lifesciences-executor:

Google Life Sciences
====================

`Google Cloud Life Sciences <https://cloud.google.com/life-sciences>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for the Life Sciences API that allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as pipelines (requires Nextflow ``20.01.0`` or later).

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover, the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor, set the property ``process.executor = 'google-lifesciences'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-machineType`
* :ref:`process-memory`
* :ref:`process-time`


See the :ref:`Google Life Sciences <google-lifesciences>` page for further configuration details.

.. _hyperqueue-executor:

HyperQueue
==========

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

The ``hyperqueue`` executor allows you to run your pipeline script by using the `HyperQueue <https://github.com/It4innovations/hyperqueue>`_ job scheduler.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``hq`` command line tool.

The pipeline must be launched from a node where the ``hq`` command is available. In a
common usage scenario, that is the cluster `head` node.

To enable the HTCondor executor, simply set ``process.executor = 'hyperqueue'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-time`


.. _htcondor-executor:

HTCondor
========

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

The ``condor`` executor allows you to run your pipeline script by using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``condor_submit`` command.

The pipeline must be launched from a node where the ``condor_submit`` command is available. In a
common usage scenario, that is the cluster `head` node.

.. note::
  The HTCondor executor for Nextflow does not currently support the HTCondor ability to transfer input/output data to
  the corresponding job computing node. Therefore, the data needs to be made accessible to the computing nodes using
  a shared file system directory from where the Nextflow workflow is executed (or specified via the ``-w`` option).

To enable the HTCondor executor, simply set ``process.executor = 'condor'`` in the ``nextflow.config`` file.

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

To enable this executor, set ``process.executor = 'ignite'`` in the ``nextflow.config`` file.

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
* :ref:`process-disk`
* :ref:`process-memory`
* :ref:`process-pod`
* :ref:`process-time`

See the :ref:`Kubernetes <k8s-page>` page to learn how to set up a Kubernetes cluster to run Nextflow pipelines.


.. _local-executor:

Local
=====

The ``local`` executor is used by default. It runs the pipeline processes on the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads`, taking advantage of the multi-core
architecture of the CPU.

The `local` executor is useful to develop and test your pipeline script on your computer, before
switching to a cluster facility when you need to run it on production data.


.. _lsf-executor:

LSF
===

The ``lsf`` executor allows you to run your pipeline script using a `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``bsub`` command.

The pipeline must be launched from a node where the ``bsub`` command is available. In a common usage
scenario, that is the cluster `head` node.

To enable the LSF executor, simply set ``process.executor = 'lsf'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. note::

    LSF supports both *per-core* and *per-job* memory limits. Nextflow assumes that LSF works in the
    *per-core* memory limits mode, thus it divides the requested :ref:`process-memory` by the number of requested :ref:`process-cpus`.

    This is not required when LSF is configured to work in the *per-job* memory limit mode. You need to specify this by
    adding the option ``perJobMemLimit`` in :ref:`config-executor` in the Nextflow configuration file.

    See also the `Platform LSF documentation <https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita>`_.


.. _moab-executor:

Moab
====

The ``moab`` executor allows you to run your pipeline script using the
`Moab <https://en.wikipedia.org/wiki/Moab_Cluster_Suite>`_ resource manager by
`Adaptive Computing <http://www.adaptivecomputing.com/>`_.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``msub`` command provided
by the resource manager.

The pipeline must be launched from a node where the ``msub`` command is available. In a common usage
scenario, that is the compute cluster `login` node.

To enable the `Moab` executor, simply set ``process.executor = 'moab'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _nqsii-executor:

NQSII
=====

The ``nsqii`` executor allows you to run your pipeline script using the `NQSII <https://www.rz.uni-kiel.de/en/our-portfolio/hiperf/nec-linux-cluster>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``qsub`` command provided
by the scheduler.

The pipeline must be launched from a node where the ``qsub`` command is available. In a common usage
scenario, that is the cluster `login` node.

To enable the NQSII executor, simply set ``process.executor = 'nqsii'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`


.. _oar-executor:

OAR
===

The ``oar`` executor allows you to run your pipeline script using the `OAR <https://oar.imag.fr>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``oarsub`` command.

The pipeline must be launched from a node where the ``oarsub`` command is available. In a common usage scenario, that is the cluster `head` node.

To enable the OAR executor, simply set ``process.executor = 'oar'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

**Known Limitations**

* Multiple ``clusterOptions`` should be semicolon-separated. This ensures that the `OAR` batch script is accurately formatted::

    clusterOptions = '-t besteffort;--project myproject'


.. _pbs-executor:

PBS/Torque
==========

The ``pbs`` executor allows you to run your pipeline script using a resource manager from the `PBS/Torque <http://en.wikipedia.org/wiki/Portable_Batch_System>`_ family of batch schedulers.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``qsub`` command provided
by the scheduler.

The pipeline must be launched from a node where the ``qsub`` command is available. In a common usage
scenario, that is the cluster `login` node.

To enable the PBS executor, simply set ``process.executor = 'pbs'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. tip::
  As of Nextflow version 23.02.0-edge or later, it is possible to specify resource settings with both the ``clusterOptions`` and
  the ``cpus`` directives by specifying the cluster options dynamically::

    cpus = 2
    clusterOptions = { "-l nodes=1:ppn=${task.cpus}:..." }

  This technique allows you to specify ``clusterOptions`` once for all processes, including any options that are specific
  to your cluster, and use the standard resource directives throughout the rest of your pipeline.


.. _pbspro-executor:

PBS Pro
=======

The ``pbspro`` executor allows you to run your pipeline script using the `PBS Pro <https://www.pbspro.org/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``qsub`` command provided
by the scheduler.

The pipeline must be launched from a node where the ``qsub`` command is available. In a common usage
scenario, that is the cluster `login` node.

To enable the PBS Pro executor, simply set ``process.executor = 'pbspro'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. tip::
  As of Nextflow version 23.02.0-edge or later, it is possible to specify resource settings with both the ``clusterOptions`` and
  the ``cpus`` and ``memory`` directives by specifying the cluster options dynamically::

    cpus = 2
    memory = 8.GB
    clusterOptions = { "-l select=1:ncpus=${task.cpus}:mem=${task.memory.toMega()}mb:..." }

  This technique allows you to specify ``clusterOptions`` once for all processes, including any options that are specific
  to your cluster, and use the standard resource directives throughout the rest of your pipeline.


.. _sge-executor:

SGE
===

The ``sge`` executor allows you to run your pipeline script using a `Sun Grid Engine <http://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_
cluster, or a compatible platform (`Open Grid Engine <http://gridscheduler.sourceforge.net/>`_, `Univa Grid Engine <http://www.univa.com/products/grid-engine.php>`_, etc).

Nextflow manages each process as a separate grid job that is submitted to the cluster using the ``qsub`` command.

The pipeline must be launched from a node where the ``qsub`` command is available. In a common usage
scenario, that is the cluster `head` node.

To enable the SGE executor, simply set ``process.executor = 'sge'`` in the ``nextflow.config`` file.

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

The ``slurm`` executor allows you to run your pipeline script using the `SLURM <https://slurm.schedmd.com/documentation.html>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster using the ``sbatch`` command.

The pipeline must be launched from a node where the ``sbatch`` command is available. In a common usage
scenario, that is the cluster `head` node.

To enable the SLURM executor, simply set ``process.executor = 'slurm'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. note:: SLURM `partitions` are comparable to job queues. Nextflow allows you to set partitions using the ``queue``
    directive listed above.

.. tip:: Nextflow does not provide direct support for SLURM multi-clusters. If you need to
  submit workflow executions to a cluster other than the current one, specify it using the
  ``SLURM_CLUSTERS`` variable in the launch environment.
