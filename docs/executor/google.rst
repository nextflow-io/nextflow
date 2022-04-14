.. _google-lifesciences-executor:

********************
Google Life Sciences
********************

Nextflow provides built-in support for Google Cloud Life Sciences API, which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions to Google Cloud.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file. Moreover the pipeline work directory must be located in a Google Storage
bucket.

To enable this executor set the property ``process.executor = 'google-lifesciences'`` in the ``nextflow.config`` file.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-disk`
* :ref:`process-machineType`
* :ref:`process-memory`

.. warning::
  This API works well for coarse-grained workloads, i.e. long-running jobs, but is not ideal for pipelines that spawn many short-lived tasks.

Read the :ref:`Google Cloud <google-page>` page to learn more about using Nextflow with Google Cloud.
