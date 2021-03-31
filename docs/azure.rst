.. _azure-page:

************
Azure Cloud
************

Requirements
============

The support for Azure Cloud requires Nextflow version ``21.02.0-edge`` or later. If you don't have it installed
use the following command to download it in your computer::

    export NXF_EDGE=1
    curl get.nextflow.io | bash
    ./nextflow -self-update

Also the support for Azure Cloud requires adding the following setting at 
the beginning of your ``nextflow.config`` file::

  plugins { 
    id 'nf-azure'
  }


.. _azure-blobstorage:

Azure Blob Storage
===================

Nextflow has built-in support for `Azure Blob Storage <https://azure.microsoft.com/en-us/services/storage/blobs/>`_.    
Files stored in a Azure blob container can be accessed transparently in your pipeline script like any other file
in the local file system.

The Blob storage account `name` and `key` needs to be provided in the Nextflow configuration file as shown below::

    azure {
      storage {
        accountName = "<YOUR BLOB ACCOUNT NAME>"
        accountKey = "<YOUR BLOB ACCOUNT KEY>"
      }
    }

As an alternative to the account key it can also used `Shared Access Token` using the setting ``sasToken`` in place
of the ``accountKey`` attribute.

.. tip::
    When creating the `Shared Access Token` make sure to allow the resource types `Container` and `Object` and allow
    the permissions: `Read`, `Write`, `Delete`, `List`, `Add`, `Create`.

.. tip::
    The value of ``sasToken`` is the token stripped by the character ``?`` from the beginning of the token.


Once the Blob Storage credentials are set you can access the files in the blob container like local files prepending
the file path with the ``az://`` prefix followed by the container name. For example, having a container named ``my-data``
container a file named ``foo.txt`` you can access it in your Nextflow script using the following fully qualified
file path ``az://my-data/foo.txt``.

.. _azure-batch:

Azure Batch
============

`Azure Batch <https://docs.microsoft.com/en-us/azure/batch/>`_ is a managed computing service that allows the execution
of containerised workloads in the Azure cloud infrastructure.

Nextflow provides a built-in support for Azure Batch which allows the seamless deployment of a Nextflow pipeline in the cloud
offloading the process executions as Batch jobs.

Get started
-------------

1 - Create a Batch account in Azure portal. Take note of the account name and key.

2 - Make sure to adjust your quotas on the pipeline's needs. There are limits on certain resources associated with the
Batch account. Many of these limits are default quotas applied by Azure at the subscription or account level.
Quotas impact on the number of Pools, CPUs and Jobs you can create at any given time.

3 - Create a Storage account and, within, an Azure Blob Container in the same location where the Batch account was created.
Take note of the account name and key.

4 - Associate the Blob Storage account with the Azure Batch account.

5 - Make sure your pipeline processes specify one or more Docker containers by using the :ref:`process-container` directive.

6 - The container images need to be published into Docker registry such as `Docker Hub <https://hub.docker.com/>`_,
`Quay <https://quay.io/>`_ or `Azure Container Registry <https://docs.microsoft.com/en-us/azure/container-registry/>`_
that can be reached by Azure Batch environment.


A minimal configuration looks like the following snippet::

    plugins {
      id 'nf-azure'
    }

    process {
      executor = 'azurebatch'
    }

    azure {
      storage {
        accountName = "<YOUR STORAGE ACCOUNT NAME>"
        accountKey = "<YOUR STORAGE ACCOUNT KEY>"
      }
      batch {
        location = '<YOUR LOCATION>'
        accountName = '<YOUR BATCH ACCOUNT NAME>'
        accountKey = '<YOUR BATCH ACCOUNT KEY>'
        autoPoolMode = true
      }
    }

In the above example, replace the location and the account placeholders with the value corresponding to your configuration and
save it to a file named ``nextflow.config``.

Given the previous configuration, launch the execution of the pipeline using the following command::

    nextflow run <PIPELINE NAME> -w az://YOUR-CONTAINER/work


Replacing ``<PIPELINE NAME>`` with a pipeline name e.g. ``nextflow-io/rnaseq-nf`` and ``YOUR-CONTAINER`` a blob
container in the storage account defined in the above configuration.

See the `Batch documentation <https://docs.microsoft.com/en-us/azure/batch/quick-create-portal>`_ for further
details about the configuration for the Azure Batch service.


Pools configuration
-------------------

When using the ``autoPoolMode`` the setting Nextflow automatically creates a `pool` of computing nodes to execute the
jobs run by your pipeline. By default it only uses 1 compute node of ``Standard_D4_v3`` type.

The pool is not removed when the pipeline execution terminates, unless the configuration setting ``deletePoolsOnCompletion=true``
is added in your pipeline configuration file.

.. warning::
    Don't forget to clean up the Batch pools to avoid in extra charges in the Batch account or use the auto scaling feature.

.. warning::
   Make sure your Batch account has enough resources to satisfy the pipeline's requirements and the pool configuration.

Pool specific settings, e.g. VM type and count, should be provided in the ``auto`` pool configuration scope, e.g. ::

    azure {
        batch {
            pools {
                auto {
                   vmType = 'Standard_D2_v2'
                   vmCount = 10
                }
            }
        }
    }


Named pools
-------------

If you want to have a more precise control on the computing nodes pools used in your pipeline using a different pool
depending on the task in your pipeline, you can use the Nextflow :ref:`process-queue` directive to specify the *name* of a
Azure Batch compute pool that has to be used to run that process' tasks.

The pool is expected to be already available in the Batch environment, unless the setting ``allowPoolCreation=true`` is
provided in the ``batch`` setting in the pipeline configuration file. In the latter case Nextflow will create the pools on-demand.

The configuration details for each pool can be specified using a snippet as shown below::

    azure {
        batch {
            pools {
                foo {
                   vmType = 'Standard_D2_v2'
                   vmCount = 10
                }

                bar {
                    vmType = 'Standard_E2_v3'
                    vmCount = 5
                }
            }
        }
    }

The above example defines the configuration for two node pools. The first will provision 10 compute nodes of type ``Standard_D2_v2``,
the second 5 nodes of type ``Standard_E2_v3``. See the `Advanced settings`_ below for the complete list of available
configuration options.

Pool autoscaling
----------------
 
Azure Batch can automatically scale pools based on parameters that you define, saving you time and money. With automatic scaling,
Batch dynamically adds nodes to a pool as task demands increase, and removes compute nodes as task demands decrease.

To enable this feature for pools created by Nextflow, add the option ``autoScale = true`` to the corresponding pool configuration scope.
For example, when using the ``autoPoolMode``, the setting looks like::

    azure {
        batch {
            pools {
                auto {
                   autoScale = true
                   vmType = 'Standard_D2_v2'
                   vmCount = 5
                   maxVmCount = 50
                }
            }
        }
    }

Nextflow uses the formula shown below to determine the number of VMs to be provisioned in the pool::

        // Get pool lifetime since creation.
        lifespan = time() - time("{{poolCreationTime}}");
        interval = TimeInterval_Minute * {{scaleInterval}};

        // Compute the target nodes based on pending tasks.
        // $PendingTasks == The sum of $ActiveTasks and $RunningTasks
        $samples = $PendingTasks.GetSamplePercent(interval);
        $tasks = $samples < 70 ? max(0, $PendingTasks.GetSample(1)) : max( $PendingTasks.GetSample(1), avg($PendingTasks.GetSample(interval)));
        $targetVMs = $tasks > 0 ? $tasks : max(0, $TargetDedicatedNodes/2);
        targetPoolSize = max(0, min($targetVMs, {{maxVmCount}}));

        // For first interval deploy 1 node, for other intervals scale up/down as per tasks.
        $TargetDedicatedNodes = lifespan < interval ? {{vmCount}} : targetPoolSize;
        $NodeDeallocationOption = taskcompletion;


The above formula initialises a pool with the number of VMs specified by the ``vmCount`` option, it scales up the pool on-demand,
based on the number of pending tasks up to ``maxVmCount`` nodes. If no jobs are submitted for execution, it scales down
to zero nodes automatically.

If you need a different strategy you can provide your own formula using the ``scaleFormula`` option.
See the `Azure Batch <https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling>`_ documentation for details.


Advanced settings
==================

The following configuration options are available:

============================================== =================
Name                                           Description
============================================== =================
azure.storage.accountName                       The blob storage account name
azure.storage.accountKey                        The blob storage account key
azure.storage.sasToken                          The blob storage shared access signature token. This can be provided as an alternative to the ``accountKey`` setting.
azure.storage.tokenDuration                     The duration of the shared access signature token created by Nextflow when the ``sasToken`` option is *not* specified (default: ``12h``).
azure.batch.accountName                         The batch service account name.
azure.batch.accountKey                          The batch service account key.
azure.batch.endpoint                            The batch service endpoint e.g. ``https://nfbatch1.westeurope.batch.azure.com``.
azure.batch.location                            The batch service location e.g. ``westeurope``. This is not needed when the endpoint is specified.
azure.batch.autoPoolMode                        Enable the automatic creation of batch pools depending on the pipeline resources demand (default: ``true``)
azure.batch.allowPoolCreation                   Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: ``false``)
azure.batch.deleteJobsOnCompletion              Enable the automatic deletion of jobs created by the pipeline execution (default: ``true``).
azure.batch.deletePoolsOnCompletion             Enable the automatic deletion of compute node pools upon pipeline completion (default: ``false``).
azure.batch.copyToolInstallMode                 Specify where the `azcopy` tool used by Nextflow. When ``node`` is specified it's copied once during the pool creation. When ``task`` is provider, it's installed for each task execution (default: ``node``)
azure.batch.pools.<name>.vmType                 Specify the virtual machine type used by the pool identified with ``<name>``.
azure.batch.pools.<name>.vmCount                Specify the number of virtual machines provisioned by the pool identified with ``<name>``.
azure.batch.pools.<name>.maxVmCount             Specify the max of virtual machine when using auto scale option.
azure.batch.pools.<name>.autoScale              Enable autoscaling feature for the pool identified with ``<name>``.
azure.batch.pools.<name>.scaleFormula           Specify the scale formula for the pool identified with ``<name>``. See Azure Batch `scaling documentation <https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling>`_ for details.
azure.batch.pools.<name>.scaleInterval          Specify the interval at which to automatically adjust the Pool size according to the autoscale formula. The minimum and maximum value are 5 minutes and 168 hours respectively (default: `10 mins`)
azure.batch.pools.<name>.schedulePolicy         Specify the scheduling policy for the pool identified with ``<name>``. It can be either ``spread`` or ``pack`` (default: ``spread``).
============================================== =================
