(azure-page)=

# Azure Cloud

:::{versionadded} 21.04.0
:::

(azure-blobstorage)=

## Azure Blob Storage

Nextflow has built-in support for [Azure Blob Storage](https://azure.microsoft.com/en-us/services/storage/blobs/). Files stored in an Azure blob container can be accessed transparently in your pipeline script like any other file in the local file system.

The Blob storage account name and key need to be provided in the Nextflow configuration file as shown below:

```groovy
azure {
    storage {
        accountName = "<YOUR BLOB ACCOUNT NAME>"
        accountKey = "<YOUR BLOB ACCOUNT KEY>"
    }
}
```

Alternatively, the **Shared Access Token** can be specified with the `sasToken` option instead of `accountKey`.

:::{tip}
When creating the Shared Access Token, make sure to allow the resource types `Container` and `Object` and allow the permissions: `Read`, `Write`, `Delete`, `List`, `Add`, `Create`.
:::

:::{tip}
The value of `sasToken` is the token stripped by the character `?` from the beginning of the token.
:::

Once the Blob Storage credentials are set, you can access the files in the blob container like local files by prepending the file path with `az://` followed by the container name. For example, a blob container named `my-data` with a file named `foo.txt` can be specified in your Nextflow script as `az://my-data/foo.txt`.

## Azure File Shares

*New in `nf-azure` version `0.11.0`*

Nextflow has built-in support also for [Azure Files](https://azure.microsoft.com/en-us/services/storage/files/). Files available in the serverless Azure File shares can be mounted concurrently on the nodes of a pool executing the pipeline. These files become immediately available in the file system and can be referred as local files within the processes. This is especially useful when a task needs to access large amounts of data (such as genome indexes) during its execution. An arbitrary number of File shares can be mounted on each pool node.

The Azure File share must exist in the storage account configured for Blob Storage. The name of the source Azure File share and mount path (the destination path where the files are mounted) must be provided. Additional mount options (see the Azure Files documentation) can be set as well for further customisation of the mounting process.

For example:

```groovy
azure {
    storage {
        accountName = "<YOUR BLOB ACCOUNT NAME>"
        accountKey = "<YOUR BLOB ACCOUNT KEY>"
        fileShares {
            <YOUR SOURCE FILE SHARE NAME> {
                mountPath = "<YOUR MOUNT DESTINATION>"
                mountOptions = "<SAME AS MOUNT COMMAND>" //optional
            }
            <YOUR SOURCE FILE SHARE NAME> {
                mountPath = "<YOUR MOUNT DESTINATION>"
                mountOptions = "<SAME AS MOUNT COMMAND>" //optional
            }
        }
    }
}
```

The files in the File share are available to the task in the directory: `<YOUR MOUNT DESTINATION>/<YOUR SOURCE FILE SHARE NAME>`.

For instance, given the following configuration:

```groovy
azure {
    storage {
        // ...

        fileShares {
            dir1 {
                mountPath = "/mnt/mydata/"
            }
        }
    }
}
```

The task can access the File share in `/mnt/mydata/dir1`.

(azure-batch)=

## Azure Batch

[Azure Batch](https://docs.microsoft.com/en-us/azure/batch/) is a managed computing service that allows the execution of containerised workloads in the Azure cloud infrastructure.

Nextflow provides built-in support for Azure Batch, allowing the seamless deployment of Nextflow pipelines in the cloud, in which tasks are offloaded as Batch jobs.

Read the {ref}`Azore Batch executor <azurebatch-executor>` section to learn more about the `azurebatch` executor in Nextflow.

### Get started

1. Create a Batch account in the Azure portal. Take note of the account name and key.
2. Make sure to adjust your quotas to the pipeline's needs. There are limits on certain resources associated with the Batch account. Many of these limits are default quotas applied by Azure at the subscription or account level. Quotas impact the number of Pools, CPUs and Jobs you can create at any given time.
3. Create a Storage account and, within that, an Azure Blob Container in the same location where the Batch account was created. Take note of the account name and key.
4. If you plan to use Azure Files, create an Azure File share within the same Storage account and upload your input data.
5. Associate the Storage account with the Azure Batch account.
6. Make sure every process in your pipeline specifies one or more Docker containers with the {ref}`process-container` directive.
7. Make sure all of your container images are published in a Docker registry that can be accessed by your Azure Batch environment, such as [Docker Hub](https://hub.docker.com/), [Quay](https://quay.io/), or [Azure Container Registry](https://docs.microsoft.com/en-us/azure/container-registry/) .

A minimal Nextflow configuration for Azure Batch looks like the following snippet:

```groovy
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
```

In the above example, replace the location placeholder with the name of your Azure region and the account placeholders with the values corresponding to your configuration.

:::{tip}
The list of Azure regions can be found by executing the following Azure CLI command:

```bash
az account list-locations -o table
```
:::

Finally, launch your pipeline with the above configuration:

```bash
nextflow run <PIPELINE NAME> -w az://YOUR-CONTAINER/work
```

Replacing `<PIPELINE NAME>` with a pipeline name e.g. `nextflow-io/rnaseq-nf` and `YOUR-CONTAINER` with a blob container in the storage account defined in the above configuration.

See the [Batch documentation](https://docs.microsoft.com/en-us/azure/batch/quick-create-portal) for further details about the configuration for Azure Batch.

### Pools configuration

When using the `autoPoolMode` option, Nextflow automatically creates a `pool` of compute nodes to execute the jobs in your pipeline. By default, it only uses one compute node of the type `Standard_D4_v3`.

The pool is not removed when the pipeline terminates, unless the configuration setting `deletePoolsOnCompletion = true` is added in your Nextflow configuration file.

Pool specific settings, such as VM type and count, should be provided in the `auto` pool configuration scope, for example:

```groovy
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
```

:::{warning}
To avoid any extra charges in the Batch account, remember to clean up the Batch pools or use auto scaling.
:::

:::{warning}
Make sure your Batch account has enough resources to satisfy the pipeline's requirements and the pool configuration.
:::

:::{warning}
Nextflow uses the same pool ID across pipeline executions, if the pool features have not changed. Therefore, when using `deletePoolsOnCompletion = true`, make sure the pool is completely removed from the Azure Batch account before re-running the pipeline. The following message is returned when the pool is still shutting down:

```
Error executing process > '<process name> (1)'
Caused by:
    Azure Batch pool '<pool name>' not in active state
```
:::

### Named pools

If you want to have more precise control over the compute node pools used in your pipeline, such as using a different pool depending on the task in your pipeline, you can use the {ref}`process-queue` directive in Nextflow to specify the ID of a Azure Batch compute pool that should be used to execute that process.

The pool is expected to be already available in the Batch environment, unless the setting `allowPoolCreation = true` is provided in the `azure.batch` config scope in the pipeline configuration file. In the latter case, Nextflow will create the pools on-demand.

The configuration details for each pool can be specified using a snippet as shown below:

```groovy
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
```

The above example defines the configuration for two node pools. The first will provision 10 compute nodes of type `Standard_D2_v2`, the second 5 nodes of type `Standard_E2_v3`. See the {ref}`Azure configuration <config-azure>` section for the complete list of available configuration options.

:::{warning}
The pool name can only contain alphanumeric, hyphen and underscore characters.
:::

:::{warning}
If the pool name includes a hyphen, make sure to wrap it with single quotes. For example::

```groovy
azure {
    batch {
        pools {
            'foo-2' {
                ...
            }
        }
    }
}
```
:::

### Requirements on pre-existing named pools

When Nextflow is configured to use a pool already available in the Batch account, the target pool must satisfy the following requirements:

1. The pool must be declared as `dockerCompatible` (`Container Type` property).
2. The task slots per node must match the number of cores for the selected VM. Otherwise, Nextflow will return an error like "Azure Batch pool 'ID' slots per node does not match the VM num cores (slots: N, cores: Y)".

### Pool autoscaling

Azure Batch can automatically scale pools based on parameters that you define, saving you time and money. With automatic scaling, Batch dynamically adds nodes to a pool as task demands increase, and removes compute nodes as task demands decrease.

To enable this feature for pools created by Nextflow, add the option `autoScale = true` to the corresponding pool configuration scope. For example, when using the `autoPoolMode`, the setting looks like:

```groovy
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
```

Nextflow uses the formula shown below to determine the number of VMs to be provisioned in the pool:

```
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
```

The above formula initialises a pool with the number of VMs specified by the `vmCount` option, and scales up the pool on-demand, based on the number of pending tasks, up to `maxVmCount` nodes. If no jobs are submitted for execution, it scales down to zero nodes automatically.

If you need a different strategy, you can provide your own formula using the `scaleFormula` option. See the [Azure Batch](https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling) documentation for details.

### Pool nodes

When Nextflow creates a pool of compute nodes, it selects:

- the virtual machine image reference to be installed on the node
- the Batch node agent SKU, a program that runs on each node and provides an interface between the node and the Batch service

Together, these settings determine the Operating System and version installed on each node.

By default, Nextflow creates pool nodes based on CentOS 8, but this behavior can be customised in the pool configuration. Below are configurations for image reference/SKU combinations to select two popular systems.

- Ubuntu 20.04 (default):

  ```groovy
  azure.batch.pools.<name>.sku = "batch.node.ubuntu 20.04"
  azure.batch.pools.<name>.offer = "ubuntu-server-container"
  azure.batch.pools.<name>.publisher = "microsoft-azure-batch"
  ```

- CentOS 8:

  ```groovy
  azure.batch.pools.<name>.sku = "batch.node.centos 8"
  azure.batch.pools.<name>.offer = "centos-container"
  azure.batch.pools.<name>.publisher = "microsoft-azure-batch"
  ```

In the above snippet, replace `<name>` with the name of your Azure node pool.

See the {ref}`Azure configuration <config-azure>` section and the [Azure Batch nodes](https://docs.microsoft.com/en-us/azure/batch/batch-linux-nodes) documentation for more details.

### Private container registry

:::{versionadded} 21.05.0-edge
:::

A private container registry for Docker images can be specified as follows:

```groovy
azure {
    registry {
        server = '<YOUR REGISTRY SERVER>' // e.g.: docker.io, quay.io, <ACCOUNT>.azurecr.io, etc.
        userName = '<YOUR REGISTRY USER NAME>'
        password = '<YOUR REGISTRY PASSWORD>'
    }
}
```

The private registry is an addition, not a replacement, to the existing configuration. Public images from other registries will still be pulled as normal, if they are requested.

:::{note}
When using containers hosted in a private registry, the registry name must also be provided in the container name specified via the {ref}`container <process-container>` directive using the format: `[server]/[your-organization]/[your-image]:[tag]`. Read more about fully qualified image names in the [Docker documentation](https://docs.docker.com/engine/reference/commandline/pull/#pull-from-a-different-registry).
:::

### Virtual Network

:::{versionadded} 23.03.0-edge
:::

Sometimes it might be useful to create a pool in an existing [Virtual Network](https://learn.microsoft.com/en-us/azure/virtual-network/). To do so, the 
`virtualNetwork` option can be added to the pool settings as follows:

```groovy
azure {
    batch {
        pools {
            auto {
                autoScale = true
                vmType = 'Standard_D2_v2'
                vmCount = 5
                virtualNetwork = '<YOUR SUBNET ID>'
            }
        }
    }
}
```

The value of the setting must be the identifier of a subnet available in the virtual network to join. A valid subnet ID has the following form:

```
/subscriptions/<YOUR SUBSCRIPTION ID>/resourceGroups/<YOUR RESOURCE GROUP NAME>/providers/Microsoft.Network/virtualNetworks/<YOUR VIRTUAL NETWORK NAME>/subnets/<YOUR SUBNET NAME>
```

:::{warning}
Batch Authentication with Shared Keys does not allow to link external resources (like Virtual Networks) to the pool. Therefore, Active Directory Authentication must be used in conjunction with the `virtualNetwork` setting.
:::

## Active Directory Authentication

:::{versionadded} 22.11.0-edge
:::

[Service Principal](https://learn.microsoft.com/en-us/azure/active-directory/develop/howto-create-service-principal-portal) credentials can optionally be used instead of Shared Keys for Azure Batch and Storage accounts.

The Service Principal should have the at least the following role assignments:

1. Contributor
2. Storage Blob Data Reader
3. Storage Blob Data Contributor

:::{note}
To assign the necessary roles to the Service Principal, refer to the [official Azure documentation](https://learn.microsoft.com/en-us/azure/role-based-access-control/role-assignments-portal?tabs=current).
:::

The credentials for Service Principal can be specified as follows:

```groovy
azure {
    activeDirectory {
        servicePrincipalId = '<YOUR SERVICE PRINCIPAL CLIENT ID>'
        servicePrincipalSecret = '<YOUR SERVICE PRINCIPAL CLIENT SECRET>'
        tenantId = '<YOUR TENANT ID>'
    }

    storage {
        accountName = '<YOUR STORAGE ACCOUNT NAME>'
    }

    batch {
        accountName = '<YOUR BATCH ACCOUNT NAME>'
        location = '<YOUR BATCH ACCOUNT LOCATION>'
    }
}
```

## Advanced configuration

Read the {ref}`Azure configuration<config-azure>` section to learn more about advanced configuration options.
