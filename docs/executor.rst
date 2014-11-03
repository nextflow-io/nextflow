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

Local executor
===============

The `local` executor is used by default. It runs the pipeline processes in the computer where Nextflow
is launched. The processes are parallelised by spawning multiple `threads` and by taking advantage of multi-cores
architecture provided by the CPU.

In a common usage scenario, the `local` executor can be useful to develop and test your pipeline script in your computer,
switching to a cluster facility when you need to run it on production data.


.. _sge-executor:

SGE executor
=============

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

LSF executor
==============

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



.. _slurm-executor:

SLURM executor
================


The `SLURM` executor allows you to run your pipeline script by using the `SLURM <https://computing.llnl.gov/linux/slurm/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``sbatch`` command.

Being so, the pipeline must be launched from a node where the ``sbatch`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the SLURM executor simply set to ``process.executor`` property to ``slurm`` value in the ``nextflow.config`` file.

The amount of resources requested by each job submission is defined by the following process directives:

* :ref:`process-cpus`
* :ref:`process-time`
* :ref:`process-memory`
* :ref:`process-clusterOptions`


.. _pbs-executor:

PBS/Torque executor
====================

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


.. _drmaa-executor:

DRMAA executor
====================

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
=========

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








