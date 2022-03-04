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
