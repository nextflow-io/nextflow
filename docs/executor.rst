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


The `SLURM` executor allows you to run your pipeline script by using the `SLURM <https://slurm.schedmd.com/documentation.html>`_ resource manager.

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

.. tip:: Nextflow does not provide a direct support for SLURM multi-clusters feature. If you need to
  submit workflow executions to a cluster that is not the current one, specify it setting the
  ``SLURM_CLUSTERS`` variable in the launching environment. 

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

.. _k8s-executor:

Kubernetes
==========

Nextflow provides an experimental support for `Kubernetes <http://kubernetes.io/>`_ clustering technology. It allows
you to deploy and transparently run a Nextflow pipeline in a Kubernetes cluster.

The following directives can be used to define the amount of computing resources needed and the container(s) to use:

* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-container`

See the :ref:`Kubernetes documentation <k8s-page>` to learn how to deploy a workflow execution in a Kubernetes cluster.

.. _awsbatch-executor:

AWS Batch
==========

Nextflow supports `AWS Batch <https://aws.amazon.com/batch/>`_ service which allows submitting jobs in the cloud
without having to spin out and manage a cluster of virtual machines. AWS Batch uses Docker containers to run tasks,
which makes deploying pipelines much simpler.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'awsbatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads. Moreover a S3 bucket must be used as pipeline work directory.

See the :ref:`AWS Batch<awscloud-batch>` page for further configuration details.

.. _google-pipelines-executor:

Google Pipelines
================

`Genomics Pipelines <https://cloud.google.com/genomics/>`_ is a managed computing service that allows the execution of
containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Genomics Pipelines API which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as pipelines (it requires Nextflow 19.1.0 or later).

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor set the property ``process.executor = 'google-pipelines'`` in the ``nextflow.config`` file.

See the :ref:`Google Pipelines <google-pipelines>` page for further configuration details.

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
    

Then you will be able to run your workflow over TES using the usual Nextflow command line, i.e.::

    nextflow run rnaseq-nf

.. note:: If the variable ``NXF_EXECUTOR_TES_ENDPOINT`` is omitted the default endpoint is ``http://localhost:8000``.

.. tip:: You can use a local `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ server using the following launch
  command line::

  ./funnel server --Server.HTTPPort 8000 --LocalStorage.AllowedDirs $HOME run

  (tested with version 0.8.0 on macOS)

.. warning:: Make sure the TES back-end can access the workflow work directory when
  data is exchanged using a local or shared file system.


**Known limitation**

* Automatic deployment of workflow scripts in the `bin` folder is not supported.
* Process output directories are not supported. For details see `#76 <https://github.com/ga4gh/task-execution-schemas/issues/76>`_.
* Glob patterns in process output declarations are not supported. For details see `#77 <https://github.com/ga4gh/task-execution-schemas/issues/77>`_.