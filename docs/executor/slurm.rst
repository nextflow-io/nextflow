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

.. note::
    SLURM partitions can be considered queues, i.e. Nextflow allows you to set partitions by using the ``queue``
    directive.

.. note::
    Nextflow does not provide direct support for SLURM's multi-clusters feature. If you need to
    submit workflow executions to a cluster that is not the current one, specify it by setting the
    ``SLURM_CLUSTERS`` variable in the launch environment.
