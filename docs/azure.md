(azure-page)=

# Azure

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

The files in the File share are available to the task in the directory: `<YOUR MOUNT DESTINATION>`.

For instance, given the following configuration:

```groovy
azure {
    storage {
        // ...

        fileShares {
          rnaseqResources {
                mountPath = "/mnt/mydata/myresources"
            }
        }
    }
}
```

The task can access the File share in `/mnt/mydata/myresources`. Note: The string `rnaseqResources` in the above config can be any name of your choice, and it does not affect the underlying mount.

:::{warning}
Azure File shares do not support authentication and management with Active Directory. The storage account key must be
set in the configuration if a share is mounted.
:::

(azure-batch)=

## Azure Batch

:::{tip}
This section describes how to manually set up and use Nextflow with Azure Batch.
You may be interested in using [Batch Forge](https://docs.seqera.io/platform/latest/compute-envs/azure-batch#compute-environment) in [Seqera Platform](https://seqera.io/platform/),
which automatically creates the required Azure infrastructure for you with minimal intervention.
:::

[Azure Batch](https://docs.microsoft.com/en-us/azure/batch/) is a managed computing service that allows the execution of containerised workloads in the Azure cloud infrastructure.

Nextflow provides built-in support for Azure Batch, allowing the seamless deployment of Nextflow pipelines in the cloud, in which tasks are offloaded as Batch jobs.

Read the {ref}`Azure Batch executor <azurebatch-executor>` section to learn more about the `azurebatch` executor in Nextflow.

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

When using the `autoPoolMode` option, Nextflow automatically creates a `pool` of compute nodes appropriate for your pipeline.

By default, the `cpus` and `memory` directives are used to find the smallest machine type that fits the requested resources in the Azure machine family, specified by `machineType`. If `memory` is not specified, 1 GB of memory is allocated per CPU. When no options are specified, it only uses one compute node of the type `Standard_D4_v3`.

To specify multiple Azure machine families, use a comma separated list with glob (`*`) values in the `machineType` directive. For example, the following will select any machine size from D or E v5 machines, with additional data disk, denoted by the `d` suffix:

```groovy
process.machineType = "Standard_D*d_v5,Standard_E*d_v5"
```

For example, the following process will create a pool of `Standard_E4d_v5` machines based when using `autoPoolMode`:

```nextflow
process EXAMPLE_PROCESS {
    machineType "Standard_E*d_v5"
    cpus 16
    memory 8.GB

    script:
    """
    echo "cpus: ${task.cpus}"
    """
}
```

Note when creating tasks that use fewer than 4 CPUs, Nextflow will create a pool with machines that have 4 times the number of CPUs required in order to pack more tasks onto each machine. This means the pipeline spends less time waiting for machines to be created, startup and join the Azure Batch pool. Similarly, if a process requires fewer than 8 CPUs Nextflow will use a machine with double the number of CPUs required. If you wish to override this behaviour you can use a specific `machineType` directive, e.g. using a `machineType` directive of `Standard_E2d_v5` will use always use a Standard_E2d_v5 machine.

The pool is not removed when the pipeline terminates, unless the configuration setting `deletePoolsOnCompletion = true` is added in your Nextflow configuration file.

Pool specific settings should be provided in the `auto` pool configuration scope. If you wish to specify a single machine size for all processes, you can specify a fixed `vmSize` for the `auto` pool.

```groovy
azure {
    batch {
        pools {
            auto {
                vmType = 'Standard_D2_v2'
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
If the pool name includes a hyphen, make sure to wrap it with single quotes. For example:

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
3. Unless you are using [Fusion](./fusion.md), all tasks must have [AzCopy](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10) available in the path. If `azure.batch.copyToolInstallMode = 'node'` this will require every node to have the azcopy binary located at `$AZ_BATCH_NODE_SHARED_DIR/bin/`.

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

### Hybrid workloads

Nextflow allows the use of multiple executors in the same workflow application. This feature enables the deployment of hybrid workloads in which some jobs are executed in the local computer or local computing cluster and some jobs are offloaded to Azure Batch.

To enable this feature, use one or more {ref}`config-process-selectors` in your Nextflow configuration to apply the Azure Batch configuration to the subset of processes that you want to offload. For example:

```groovy
process {
    withLabel: bigTask {
        executor = 'azurebatch'
        queue = 'my-batch-pool'
        container = 'my/image:tag'
    }
}

azure {
    storage {
        accountName = '<YOUR STORAGE ACCOUNT NAME>'
        accountKey = '<YOUR STORAGE ACCOUNT KEY>'
    }
    batch {
        location = '<YOUR LOCATION>'
        accountName = '<YOUR BATCH ACCOUNT NAME>'
        accountKey = '<YOUR BATCH ACCOUNT KEY>'
    }
}
```

With the above configuration, processes with the bigTask {ref}`process-label` will run on Azure Batch, while the remaining processes will run on the local computer.

Then launch the pipeline with the `-bucket-dir` option to specify an Azure Blob Storage path for the jobs computed with Azure Batch, and optionally, use the `-work-dir` option to specify the local storage for the jobs computed locally:

```bash
nextflow run <script or project name> -bucket-dir az://my-container/some/path
```

:::{warning}
The Azure Blob Storage path needs to contain at least one sub-directory (e.g. `az://my-container/work` rather than `az://my-container`).
:::

:::{note}
Nextflow will automatically manage the transfer of input and output files between the local and cloud environments when using hybrid workloads.

:::{tip}
When using [Fusion](./fusion.md), the `-bucket-dir` option is not required. Fusion implements a distributed virtual file system that allows seamless access to Azure Blob Storage using a standard POSIX interface, enabling direct mounting of remote blob storage as if it were a local file system. This simplifies and speeds up most operations, bridging the gap between cloud-native storage and data analysis workflows.
:::

## Microsoft Entra

Using Microsoft Entra for role-based access control is more secure than using access keys and should be used wherever possible. You can authenticate to Azure Entra using a Managed Identity when running on resources within the Azure environment, or by authenticating as an Azure Service Principal when running on external resources.

### Required role assignments

To access Azure resources, you must have the relevant role assignments (permissions). Access Azure Blob storage data (e.g. to retrieve data from a private Azure Storage account) you must have the following permissions:

1. Storage Blob Data Reader
2. Storage Blob Data Contributor

To run Nextflow on Azure Batch you must have the following permissions to create and destroy resources in the Azure Batch account:

1. Batch Contributor

To assign the necessary roles to a Managed Identity or Service Principal, refer to the [official Azure documentation](https://learn.microsoft.com/en-us/azure/role-based-access-control/role-assignments-portal?tabs=current).

(azure-managed-identities)=

### Managed identities

:::{versionadded} 24.05.0-edge
:::

An Azure [Managed Identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) can be used to authenticate with Azure Resources without the use of access keys or credentials. When using a Managed Identity, an Azure resource is able to authenticate because of what _it is_, instead of using access keys. For example, if Nextflow is running on an [Azure Virtual Machine with a managed identity enabled](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-to-configure-managed-identities?pivots=qs-configure-portal-windows-vm) which had the relevant permissions, it would be able to run a Nextflow workflow on Azure Batch and use data from a private storage account with no additional credentials supplied. This is more secure and reliable than using Access Keys or a Service Principal which can be compromised. The limitation is they are only able to work from within the Azure account, i.e. you cannot use a managed identity from an external service.

An Azure Managed identity comes in two forms, system-assigned and user-assigned. Both are functionally equivalent but have a slightly method

#### System Assigned Managed Identity

When running on an Azure Service such as an Azure Virtual Machine, you can enable system-assigned Managed Identity for that machine. This grants the machine an identity in Azure from which it receives the permissions and role assignments.

1. First we must [enable system-assigned Managed Identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-to-configure-managed-identities?pivots=qs-configure-portal-windows-vm). Once you have done this the machine has an identity in Azure Entra.
2. Next we must add the relevant role assignments to this Managed Identity. On the Azure Portal page for the virtual machine, select 'Identity' and then click 'Azure Role Assignments' to modify the existing role assignments. Note you must have `Microsoft.Authorization/roleAssignments/write` to perform this action.
3. Make sure the identity has the following role assignments:
    - Storage Blob Data Reader
    - Storage Blob Data Contributor
    - Batch Contributor
4. Save the changes.
5. Use the following Nextflow configuration to enable Nextflow to adopt the system-assigned identity while running on this machine:

```groovy
process.executor = 'azurebatch'
azure {
    managedIdentity {
        system = true
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

#### User Assigned Managed Identity

A system-assigned managed identity is essentially 'anonymous' and is tied to a single resource. By comparison, a user-assigned managed identity is created by the user and can be assigned to multiple resources, furthermore the lifecycle of a user-assigned managed identity is not tied to the resource. See [the Azure Documentation](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/managed-identity-best-practice-recommendations#choosing-system-or-user-assigned-managed-identities) for further details.

We can add a user-assigned identity to a resource in a similar manner to a system-assigned identity, but we must create it first.

1. Create a Managed Identity as per [the Azure Documentation](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-manage-user-assigned-managed-identities?pivots=identity-mi-methods-azp).
2. Assign the relevant role permissions to the Managed Identity as before.
3. Assign the Managed Identity to the Azure Resource, this can be done at creation time or afterwards. [As an example, see the documentation on how to do this for a virtual machine](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-to-configure-managed-identities?pivots=qs-configure-portal-windows-vm#user-assigned-managed-identity).
4. Retrieve the client ID for the Managed Identity. On the Azure Portal, this can be found on the 'Overview' or 'Properties' page as 'client ID'.
5. Use the following configuration to enable Nextflow to adopt the user-assigned identity while running on the Azure Resource:

```groovy
process.executor = 'azurebatch'
azure {
    managedIdentity {
        clientId = '<USER ASSIGNED MANAGED IDENTITY CLIENT ID>'
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

(azure-service-principal)=

### Service Principals

:::{versionadded} 22.11.0-edge
:::

[Service Principal](https://learn.microsoft.com/en-us/azure/active-directory/develop/howto-create-service-principal-portal) credentials can be used access to Azure Batch and Storage accounts. Similar to a Managed Identity, a Service Principal is an account which can have specific permissions and role based access. However, unlike with Managed Identities you must use a secret key to authenticate as a Service Principal. However, this means you can access authenticate as a Service Principal when operating outside of the Azure account.

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
