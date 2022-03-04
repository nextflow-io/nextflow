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
