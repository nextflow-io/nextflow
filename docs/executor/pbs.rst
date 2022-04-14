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
