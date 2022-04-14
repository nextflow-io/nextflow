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
