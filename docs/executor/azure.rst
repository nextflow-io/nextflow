.. _azurebatch-executor:

***********
Azure Batch
***********

Nextflow provides built-in support for Azure Batch, which allows the seamless deployment of a Nextflow pipeline in the cloud,
offloading the process executions as Batch jobs.

The pipeline processes must specify the Docker image to use by defining the ``container`` directive, either in the pipeline
script or the ``nextflow.config`` file.

To enable this executor set the property ``process.executor = 'azurebatch'`` in the ``nextflow.config`` file.

The pipeline can be launched either in a local computer or a cloud virtual machine. The latter is suggested for heavy or long
running workloads. Moreover a Azure Blob storage container must be used as pipeline work directory.

Read the :ref:`Azore Cloud <azure-page>` page to learn more about using Nextflow with Azure Cloud.
