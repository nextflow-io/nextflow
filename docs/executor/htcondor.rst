.. _htcondor-executor:

HTCondor
========

.. warning:: This feature is experimental, and may change in a future release.

The ``condor`` executor allows you to run your pipeline script by using the `HTCondor <https://research.cs.wisc.edu/htcondor/>`_ resource manager.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``condor_submit`` command.

Being so, the pipeline must be launched from a node where the ``condor_submit`` command is available, that is, in a
common usage scenario, the cluster `head` node.

.. note::
  The HTCondor executor for Nextflow does not currently support HTCondor's ability to transfer input/output data to
  the corresponding job compute node. Therefore the data needs to be made accessible to the compute nodes using
  a shared file system directory from which the Nextflow workflow must be executed (or specified via the ``-w`` option).

To enable the HTCondor executor simply set ``process.executor = 'condor'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-memory`
* :ref:`process-time`
