.. _awsbatch-executor:

*********
AWS Batch
*********

Nextflow provides built-in support for AWS Batch, which allows the seamless deployment of a Nextflow pipeline
in the cloud, offloading the process executions as Batch jobs.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'awsbatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a EC2 instance. The latter is suggested for heavy or long
running workloads. Moreover, a S3 bucket must be used as the pipeline work directory.

Resource requests and other job characteristics can be controlled via the following process directives:

* :ref:`process-accelerator`
* :ref:`process-cpus`
* :ref:`process-memory`
* :ref:`process-queue`

Read the :ref:`AWS<aws-page>` page to learn more about using Nextflow with AWS.
