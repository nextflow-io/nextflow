.. _bridge-executor:

Bridge
======

`Bridge <https://github.com/cea-hpc/bridge>`_ is an abstraction layer to ease batch system and resource manager usage in 
heterogeneous HPC environments.

It is open source software and can be installed on top of existing classical job schedulers such as Slurm or LSF, or other 
schedulers. Bridge allows to submit jobs, get information on running jobs, stop jobs, get information on the cluster system, etc.

For more details on how to install the Bridge system, see the `documentation <https://github.com/cea-hpc/bridge>`_.

To enable the Bridge executor simply set ``process.executor = 'bridge'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`
