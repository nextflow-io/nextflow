.. _lsf-executor:

LSF
===

The ``lsf`` executor allows you to run your pipeline script by using a `Platform LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster.

Nextflow manages each process as a separate job that is submitted to the cluster by using the ``bsub`` command.

Being so, the pipeline must be launched from a node where the ``bsub`` command is available, that is, in a common usage
scenario, the cluster `head` node.

To enable the LSF executor simply set ``process.executor = 'lsf'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-clusterOptions`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`
* :ref:`process-time`

.. note::
    LSF supports both *per-core* and *per-job* memory limits. Nextflow assumes that LSF works in the
    *per-core* memory limit mode, thus it divides the requested :ref:`process-memory` by the number of requested :ref:`process-cpus`.

    This is not required when LSF is configured to work in *per-job* memory limit mode. You will need to specify it by
    adding the option ``perJobMemLimit`` in :ref:`config-executor` in the Nextflow configuration file.

    See also the `Platform LSF documentation <https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsf.conf.lsb_job_memlimit.5.dita>`_.
