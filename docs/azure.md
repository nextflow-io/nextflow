(azure-page)=

# Azure

## Overview

Nextflow provides built-in support for Azure Cloud Services. It enables you to:

- Store and access data using Azure Blob Storage and Azure file shares.
- Execute workflows using Azure Batch.

:::{tip}
For automated Azure infrastructure setup, consider using [Batch Forge](https://docs.seqera.io/platform/latest/compute-envs/azure-batch#compute-environment) in [Seqera Platform](https://seqera.io/platform/).
:::

## Quick start

To run pipelines with Azure Batch:

1. Create an Azure Batch account in the Azure portal.

2. Increase the quotas in your Azure Batch account to the pipeline's needs. Quotas impact the number of pools, CPUs, and jobs you can create.

3. Create a storage account and Azure Blob Storage in the same region as the Batch account.

4. Add the following settings to your Nextflow configuration file:

    - Add authentication and account details. See [Authentication](#authentication) for configuration examples.

    - Set `process.executor` to `azurebatch` to make Nextflow submit tasks to Azure Batch.

    - Set `workDir` to a working directory on Azure Blob Storage. For example, `az://<BLOB_STORAGE>/work`, where `BLOB_CONTAINER` is a blob container in your storage account.

5. Launch your pipeline with the above configuration:

    ```bash
    nextflow run <PIPELINE_NAME>
    ```

:::{tip}
You can list Azure regions with:

```bash
az account list-locations -o table
```
:::

## Authentication

Nextflow supports three authentication methods for Azure services, listed in order of security and recommended usage:

- [Managed identities](#managed-identities): The most secure option, available only within Azure. Uses Azure-managed credentials without storing secrets.

- [Service principals](#service-principals): A secure and flexible method that works across environments. Uses Microsoft Entra credentials for authentication.

- [Access keys](#access-keys): The most basic and least secure method. Relies on direct access keys for authentication.

### Required roles

The following role assignments are required to use Azure services:

For Azure Storage:

- Storage Blob Data Reader
- Storage Blob Data Contributor

For Azure Batch:

- Batch Data Contributor

To assign roles to a managed identity or service principal, see the [Azure documentation](https://learn.microsoft.com/en-us/azure/role-based-access-control/role-assignments-portal?tabs=current).

(azure-managed-identities)=

### Managed identities

:::{versionadded} 24.05.0-edge
:::

[Managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) is the recommended authentication method when running Nextflow within Azure. It automatically manages credentials without requiring you to store secrets.

When using a managed identity, an Azure resource can authenticate based on what it is rather than using access keys. For example, Nextflow can submit jobs to Azure Batch and access private storage accounts without additional credentials if it is running on Azure Virtual Machines with a managed identity and the necessary permissions.

There are two types of managed identities:

- System-assigned managed identity
- User-assigned managed identity

**System-assigned managed identity**

A system-assigned identity is tied to a specific Azure resource. To use it:

1. Enable [system-assigned managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-to-configure-managed-identities?pivots=qs-configure-portal-windows-vm#system-assigned-managed-identity) on your Azure resource.

2. Configure the required role assignments. See [Required roles](#required-roles) for more information.

3. Add the following configuration:

```groovy
azure {
    managedIdentity {
        system = true
    }
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
    }
    batch {
        accountName = '<BATCH_ACCOUNT_NAME>'
        location = '<BATCH_ACCOUNT_LOCATION>'
    }
}
```

Replace the following:

- `STORAGE_ACCOUNT_NAME`: your Azure storage account name
- `BATCH_ACCOUNT_NAME`: your Azure Batch account name
- `BATCH_ACCOUNT_LOCATION`: your Azure Batch account location

:::{tip}
Nextflow uses the `AZURE_MANAGED_IDENTITY_SYSTEM` environment variable if the managed identity is not set in the Nextflow configuration file. Set this to `true` to enable a system-assigned managed identity
:::

**User-assigned managed identity**

A user-assigned identity can be shared across multiple Azure resources and its lifecycle is not tied to any specific resource. See the [Azure documentation](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/managed-identity-best-practice-recommendations#choosing-system-or-user-assigned-managed-identities) for more information.

To use user-assigned managed identity:

1. Create a [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-manage-user-assigned-managed-identities?pivots=identity-mi-methods-azp#create-a-user-assigned-managed-identity).

2. Configure the required role assignments.

3. [Assign the identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/how-to-configure-managed-identities?pivots=qs-configure-portal-windows-vm#user-assigned-managed-identity) to your Azure resource.

4. Add the following configuration:

```groovy
azure {
    managedIdentity {
        clientId = '<USER_ASSIGNED_MANAGED_IDENTITY_CLIENT_ID>'
    }
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
    }
    batch {
        accountName = '<BATCH_ACCOUNT_NAME>'
        location = '<BATCH_ACCOUNT_LOCATION>'
    }
}
```

Replace the following:

- `USER_ASSIGNED_MANAGED_IDENTITY_CLIENT_ID`: your user assigned managed identity object ID
- `STORAGE_ACCOUNT_NAME`: your Azure Storage account name
- `BATCH_ACCOUNT_NAME`: your Azure Batch account name
- `BATCH_ACCOUNT_LOCATION`: your Azure Batch account location

:::{tip}
Nextflow uses the `AZURE_MANAGED_IDENTITY_USER` environment variable if the managed identity client ID is not provided in the Nextflow configuration file.
:::

### Service principals

:::{versionadded} 22.11.0-edge
:::

[Service principal](https://learn.microsoft.com/en-us/entra/identity-platform/app-objects-and-service-principals) credentials can be used to access Azure Batch and Storage accounts.

Similar to a managed identity, a service principal is an application identity that can have specific permissions and role-based access. Unlike managed identities, service principals require a secret key for authentication. The secret key is less secure, but it allows you to authenticate with Azure resources when operating outside of the Azure environment, such as from your local machine or a non-Azure cloud.

To create a service principal, follow the steps in the [Register a Microsoft Entra app and create a service principal](https://learn.microsoft.com/en-us/entra/identity-platform/howto-create-service-principal-portal) guide and add the credentials to your Nextflow configuration file. For example:

```groovy
azure {
    activeDirectory {
        servicePrincipalId = '<SERVICE_PRINCIPAL_CLIENT_ID>'
        servicePrincipalSecret = '<SERVICE_PRINCIPAL_CLIENT_SECRET>'
        tenantId = '<TENANT_ID>'
    }
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
    }
    batch {
        accountName = '<BATCH_ACCOUNT_NAME>'
        location = '<BATCH_ACCOUNT_LOCATION>'
    }
}
```

:::{tip}
Nextflow uses the following environment variables if storage settings are not set in the Nextflow config file:

- `AZURE_CLIENT_ID`: your Azure service principal application ID
- `AZURE_CLIENT_SECRET`: your Azure service principal secret
- `AZURE_TENANT_ID`: your Azure Entra account tenant ID
:::

### Access keys

The basic authentication method uses direct access keys. While easy to set up, it's less secure than managed identities or service principals because it provides full access to your resources. This approach requires you to obtain and manage the access keys for your Azure Batch and Storage accounts.

You can find your access keys in the Azure portal:

- For Storage accounts, navigate to your Storage account and select **Access keys** in the left-hand navigation menu.
- For Batch accounts, navigate to your Batch account and select **Keys** in the left-hand navigation menu.

Keep access keys secure and rotate them regularly. Use environment variables instead of hardcoding keys in configuration files, especially in shared or version-controlled environments:

```groovy
azure {
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
        accountKey = '<STORAGE_ACCOUNT_KEY>'
    }
    batch {
        accountName = '<BATCH_ACCOUNT_NAME>'
        accountKey = '<BATCH_ACCOUNT_KEY>'
        location = '<LOCATION>'
    }
}
```

You can also use a Shared Access Token (SAS) instead of an account key to provide time-limited access to your storage account:

```groovy
azure.storage.sasToken = '<SAS_TOKEN>'
```

:::{tip}
When creating a SAS token, make sure to enable `Read`, `Write`, `Delete`, `List`, `Add`, `Create` permissions for the `Container` and `Object` resource types. The value of `sasToken` should be stripped of the leading `?` character.
:::

:::{tip}
Nextflow uses the following environment variables if storage settings are not provided in the Nextflow configuration file:

- `AZURE_STORAGE_ACCOUNT_NAME`: your Azure Storage account name
- `AZURE_STORAGE_ACCOUNT_KEY`: your Azure Storage account key
- `AZURE_BATCH_ACCOUNT_NAME`: your Azure Batch account name
- `AZURE_BATCH_ACCOUNT_KEY`: your Azure Batch account key
- `AZURE_STORAGE_SAS_TOKEN`: your Azure Storage account SAS token
:::

## Azure Blob Storage

Azure Storage provides several services for storing and accessing data in Azure. Nextflow primarily supports Azure Blob Storage for most workloads. Limited support is available for Azure file shares for specific use cases.

Nextflow uses the `az://` URI scheme to enable transparent access to files stored in Azure Blob Storage. Once Azure Storage is configured, you can use this scheme to reference files in any blob container as if they were local. For example, the path `/path/to/file.txt` in a container named `my-container` can be accessed as follows:

```groovy
file('az://my-container/path/to/file.txt')
```

The Blob storage account name and key need to be provided in the Nextflow configuration file as shown below:

```groovy
azure {
    storage {
        accountName = '<BLOB_ACCOUNT_NAME>'
    }
}
```

:::{tip}
Nextflow uses the following environment variables if storage settings are not provided in the Nextflow configuration file:

- `AZURE_STORAGE_ACCOUNT_NAME`: your Azure Storage account name
- `AZURE_STORAGE_ACCOUNT_KEY`: your Azure Storage account access key
- `AZURE_STORAGE_SAS_TOKEN`: a shared access signature (SAS) token for Azure Storage access
:::

(azure-batch)=

## Azure Batch

[Azure Batch](https://learn.microsoft.com/en-us/azure/batch/) is a managed computing service that enables large-scale parallel and high-performance computing (HPC) batch jobs in Azure. Nextflow provides a built-in executor for Azure Batch, allowing you to offload workflow tasks to the cloud.

Nextflow integrates seamlessly with Azure Batch to:

- Dynamically create and manage compute pools based on workflow requirements.
- Automatically scale nodes up and down as needed during workflow execution.
- Support both regular and low-priority VMs for cost optimization.
- Manage task dependencies and data transfers between Azure Storage and compute nodes.

This section describes how to configure and use Azure Batch with Nextflow for efficient cloud-based workflow execution. For comprehensive information about Azure Batch features and capabilities, refer to the [official Azure Batch documentation](https://learn.microsoft.com/en-us/azure/batch/).

### Overview

Nextflow integrates with Azure Batch by mapping its execution model to Azure Batch's structure. A Nextflow process corresponds to an Azure Batch job, and every execution of that process (a Nextflow task) becomes an Azure Batch task. These Azure Batch tasks are executed on compute nodes within an Azure Batch pool, which is a collection of virtual machines that can scale up or down based on an autoscale formula.

Nextflow manages these pools dynamically. You can assign processes to specific, pre-existing pools using the process `queue` directive. Nextflow will create if it doesn't exist and `azure.batch.allowPoolCreation` is set to `true`. Alternatively, `autoPoolMode` enables Nextflow to automatically create multiple pools based on the CPU and memory requirements defined in your processes.

An Azure Batch task is created for each Nextflow task. This task first downloads the necessary input files from Azure Blob Storage to its assigned compute node. It then runs the process script. Finally, it uploads any output files back to Azure Blob Storage.

Nextflow ensures that the corresponding Azure Batch jobs are marked as terminated when the pipeline completes. The compute pools, if configured for auto-scaling, automatically scale down according to their defined formula to minimize costs. Nextflow can also be configured to delete these node pools after pipeline completion.

To use Azure Batch, set the `executor` directive to `azurebatch` in your Nextflow configuration file and add a work directory on Azure Blob Storage:

```groovy
process {
    executor = 'azurebatch'
}

workDir = 'az://<BLOB_CONTAINER>/work'

azure {
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
    }
    batch {
        location = '<LOCATION>'
        accountName = '<BATCH_ACCOUNT_NAME>'
        autoPoolMode = true
        allowPoolCreation = true
    }
}
```

:::{note}
The work directory must be a subdirectory, for example, `az://container/work`, not `az://container`.
:::

Finally, launch your pipeline with the above configuration:

```bash
nextflow run <PIPELINE_NAME>
```

### Quotas

Azure Batch enforces quotas on resources like pools, jobs, and VMs.

- **Pools:** There is a limit on the number of pools that can exist at one time. Nextflow throws an error if it attempts to create a pool when this quota is exhausted. To proceed, you must delete existing pools.

- **Active jobs:** There is a limit on the number of active jobs. Completed jobs should be terminated. If this limit is reached, Nextflow throws an error. You must terminate or delete jobs before you run additional pipelines.

- **VM cores:** Each Azure VM family has a maximum core quota. If you request a VM family size for which you lack sufficient core quota, any jobs and tasks created by Nextflow will remain pending, as the pool cannot scale due to the quota. This can cause your pipeline to hang indefinitely.

You can increase your quota on the Azure Portal by opening your Batch account and selecting **Quotas** in the left-hand menu, then request an increase in quota numbers.

### Pools

Nextflow supports two approaches for managing Batch pools: auto pools and named pools. These approaches are described in the following sections.

:::{warning}
Clean up the Batch pools or use auto-scaling to avoid additional charges in the Batch account.
:::

### Auto pools

When using the `autoPoolMode` mode, Nextflow automatically creates pools of compute nodes appropriate for your pipeline.

The `cpus` and `memory` directives are used to find the smallest machine type in the Azure machine family, specified by the `machineType` directive, that satisfies the requested resources. If memory is not specified, 1 GB of memory is requested per CPU. When no options are specified, it uses one compute node of type *Standard_D4_v3*.

To use `autoPoolMode`, enable it in the `azure.batch` config scope:

```groovy
azure {
    batch {
        autoPoolMode = true
        allowPoolCreation = true
    }
}

process {
    machineType = 'Standard_D4_v3'
    cpus = 4
    memory = 16.GB
}
```

You can specify multiple machine families per process using glob patterns in the `machineType` directive:

```groovy
// D or E v5 machines with data disk
process.machineType = "Standard_D*d_v5,Standard_E*d_v5" 
```

For example, the following process creates a pool of *Standard_E8d_v5* machines when using `autoPoolMode`:

```groovy
process EXAMPLE_PROCESS {
    machineType "Standard_E*d_v5"
    cpus 8
    memory 8.GB

    script:
    """
    echo "cpus: ${task.cpus}"
    """
}
```

Additional notes:

- For tasks using fewer than 4 CPUs, Nextflow creates pools with 4x CPUs to enable task packing.
- For tasks using fewer than 8 CPUs, Nextflow uses 2x CPUs.
- This behavior can be overridden by specifying an exact machine type (e.g., *Standard_E2d_v5*).
- Regexes can be used to avoid certain machine types (e.g., `Standard_*[^p]_v*` avoids ARM-based machines).
- The pool is not removed when the pipeline terminates unless `azure.batch.deletePoolsOnCompletion` is enabled in your Nextflow configuration file.

### Named pools

The `queue` directive can be used to specify the ID of the Azure Batch compute pool that should be used for a process.

When `azure.batch.allowPoolCreation = true` is set in your configuration, Nextflow creates the pools on-demand. Otherwise, the pool must be available in the Batch environment when Nextflow submits tasks.

Each pool can be configured separately in the Nextflow configuration. For example:

```groovy
azure {
    batch {
        pools {
            small {
                vmType = 'Standard_D2_v2'
                vmCount = 5
            }
            large {
                vmType = 'Standard_D8_v3'
                vmCount = 2
            }
        }
    }
}

process {
    withLabel: small_jobs {
        queue = 'small'
    }
    withLabel: large_jobs {
        queue = 'large'
    }
}
```

The above example defines the configuration for two node pools, `small` and `large`. The `small` pool provisions 10 compute nodes of type *Standard_D2_v2*, and the `large` pool provisions 5 nodes of type *Standard_E2_v3*. See the {ref}`config-azure` configuration scope for the list of available configuration options.

:::{warning}
Pool names may contain alphanumeric characters, hyphens, and underscores. Hyphenated names must be quoted, e.g., `'pool-1'`.
:::

**Requirements for pre-existing named pools**

The target pool must satisfy the following requirements when Nextflow is configured to use a pool that is already available in the Batch account:

- The pool must be declared as `dockerCompatible` (Container Type property).

- The task slots per node must match the number of cores for the selected VM. Otherwise, Nextflow returns an error of the form "Azure Batch pool 'ID' slots per node does not match the VM num cores (slots: N, cores: Y)".

- Unless you are using Fusion, all tasks must have the `azcopy` command line tool available in the path. If `azure.batch.copyToolInstallMode = 'node'`, every node must have `azcopy` available at `$AZ_BATCH_NODE_SHARED_DIR/bin/`.

### Auto-scaling

Azure Batch can automatically scale pools based on the parameters you define. With automatic scaling, Batch dynamically adds nodes to a pool as task demands increase and removes nodes as task demands decrease.

To enable this feature for pools created by Nextflow, add `autoScale = true` to the corresponding pool configuration scope. For example, when using the `autoPoolMode`:

```groovy
azure {
    batch {
        pools {
            <POOL_NAME> {
                autoScale = true
                vmCount = 5        // Initial count
                maxVmCount = 50    // Maximum count
            }
        }
    }
}
```

The default auto-scaling formula initializes the pool with the number of VMs specified by the `vmCount` option. It dynamically scales up the pool based on the number of pending tasks. The pool limit is defined by `maxVmCount`. The formula automatically scales the pool down to zero nodes to minimize costs when no jobs are submitted for execution. The default formula is shown below:

```
// Get pool lifetime since creation.
lifespan = time() - time("{{poolCreationTime}}");
interval = TimeInterval_Minute * {{scaleInterval}};

// Compute the target nodes based on pending tasks.
$samples = $PendingTasks.GetSamplePercent(interval);
$tasks = $samples < 70 ? max(0, $PendingTasks.GetSample(1)) : max( $PendingTasks.GetSample(1), avg($PendingTasks.GetSample(interval)));
$targetVMs = $tasks > 0 ? $tasks : max(0, $TargetDedicatedNodes/2);
targetPoolSize = max(0, min($targetVMs, {{maxVmCount}}));

// For first interval deploy 1 node, for other intervals scale up/down as per tasks.
$TargetDedicatedNodes = lifespan < interval ? {{vmCount}} : targetPoolSize;
$NodeDeallocationOption = taskcompletion;
```

The `azure.batch.pools.<POOL_NAME>.scaleFormula` setting can be used to specify custom formulas.

### Task authentication

By default, Nextflow creates SAS tokens for specific containers and passes them to tasks to enable file operations with Azure Storage. SAS tokens expire after a set period of time. The expiration time is 48 hours by default and cat be configured using `azure.storage.tokenDuration` in your configuration.

:::{versionadded} 25.05.0-edge
:::

You can also authenticate to Azure Storage using a managed identity when using Fusion.

To do this:

1. Create a user-assigned managed identity with the Azure Storage Blob Data Contributor role for your storage account.
2. Attach the user-assigned managed identity to the node pool manually.
3. Set `azure.managedIdentity.clientId` to this identity in your configuration.

```groovy
azure {
    managedIdentity {
        clientId = '<MANAGED_IDENTITY_CLIENT_ID>'
    }
}
```

Each task authenticates as the managed identity, and downloads and uploads files to Azure Storage using these credentials. It is possible to attach more than one managed identity to a pool. Fusion uses the identity specified by `azure.managedIdentity.clientId`.

### Task packing

Each Azure Batch node is allocated a specific number of task slots that determine how many tasks can run concurrently on that node. The number of slots is calculated based on the virtual machine's available CPUs. Nextflow intelligently assigns task slots to each process based on the proportion of total node resources requested through the `cpus`, `memory`, and `disk` directives.

For example, consider a *Standard_D4d_v5* machine with 4 vCPUs, 16 GB of memory, and 150 GB local disk:

- If a process requests `cpus 2`, `memory 8.GB`, or `disk 75.GB`, two task slots are allocated (50% of resources), allowing two tasks to run concurrently on the node.

- If a process requests `cpus 4`, `memory 16.GB`, or `disk 150.GB`, four task slots are allocated (100% of resources), allowing one task to run on the node.

Resource overprovisioning can occur if tasks consume more than their allocated share of resources. For instance, the node described above my become overloaded and fail if a task with `cpus 2` uses more than 8 GB of memory or 75 GB of disk space. Make sure to accurately specify resource requirements to ensure optimal performance and prevent task failures.

:::{warning}
Azure virtual machines come with fixed storage disks that are not expandable. Tasks will fail if the tasks running concurrently on a node use more storage than the machine has available.
:::

### Task containers

Every process in your pipeline must specify a container in order to be executed on Azure Batch.

Nextflow supports both public and private container images. You can pull public images from repositories such as https://community.wave.seqera.io/. You must authenticate for any private images that you use. Ensure the Batch pool has access when using a private registry like Azure Container Registry (ACR). Provide ACR credentials in the `azure.registry` configuration scope:

```groovy
azure {
    registry {
        server = '<REGISTRY_SERVER>'     // e.g., 'myregistry.azurecr.io'
        userName = '<REGISTRY_USER>'
        password = '<REGISTRY_PASSWORD>'
    }
}
```

:::{tip}
Nextflow uses the following environment variables if the registry credentials are not provided in the Nextflow configuration file:

- `AZURE_REGISTRY_USER_NAME`: the username for Azure Container Registry authentication
- `AZURE_REGISTRY_PASSWORD`: the password for Azure Container Registry authentication
:::

Container images from public registries such as Docker Hub can be used without additional configuration, even when a private registry is specified in the `azure.registry` scope. The Docker runtime on Azure Batch VMs automatically uses the appropriate registry based on the container image specified for the task.

### VM images

When Nextflow creates a pool of compute nodes, it selects:

- The virtual machine image reference (defines the operating system (OS) and software to be installed on the node).
- The Batch node agent SKU (a program that runs on each node and provides the interface between the node and the Azure Batch service).

These settings determine the operating system, version, and available features for each compute node in your pool.

By default, Nextflow creates pool nodes suitable for most bioinformatics workloads based on Ubuntu 22.04. You can customize the OS in your pool configuration to meet specific requirements.

Below are example configurations for common operating systems used in scientific computing:

**Ubuntu 22.04 (default)**

```groovy
azure {
    batch {
        pools {
            <POOL_NAME> {
                sku = "batch.node.ubuntu 22.04"
                offer = "ubuntu-hpc"
                publisher = "microsoft-dsvm"
            }
        }
    }
}
```

**CentOS 8**

```groovy
azure {
    batch {
        pools {
            <POOL_NAME> {
                sku = "batch.node.centos 8"
                offer = "centos-container"
                publisher = "microsoft-azure-batch"
            }
        }
    }
}
```

### Advanced features

**Virtual networks**

Pools can be configured to use virtual networks to connect to your existing network infrastructure.

```groovy
azure.batch.pools.<POOL_NAME>.virtualNetwork = '<SUBNET_ID>'
```

The subnet ID must be in the following format:

```
/subscriptions/<SUBSCRIPTION>/resourceGroups/<GROUP>/providers/Microsoft.Network/virtualNetworks/<VNET>/subnets/<SUBNET>
```

Nextflow can submit tasks to existing pools that are already attached to a virtual network without requiring additional configuration.

:::{warning}
Virtual networks require Microsoft Entra authentication (service principal or managed identity).
:::

**Start tasks**

Start tasks are optional commands that run when a node joins the pool. They are useful for setting up the node environment or for running other commands that need to be executed when the node is created.

```groovy
azure {
    batch {
        pools {
            <POOL_NAME> {
                startTask {
                    script = 'echo "Hello, world!"'
                    privileged = true  // optional, defaults to false
                }
            }
        }
    }
}
```

**Azure file shares**

Azure file shares provide fully managed file shares that can be mounted to compute nodes. Files become immediately available in the file system and can be accessed as local files within processes.

The Azure file share must exist in the storage account configured for Blob Storage. You must provide the name of the source Azure file share and mount path (the path where the files are mounted). Additional mount options can also be set for further customization of the mounting process. See the [Azure Files documentation](https://learn.microsoft.com/en-us/azure/storage/files/storage-files-introduction) for more information.

For example:

```groovy
azure {
    storage {
        fileShares {
            references {  // any name of your choice
                mountPath = '/mnt/references'
                mountOptions = '-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp'  // optional
            }
        }
    }
}
```

:::{warning}
File shares must be authenticated using a storage account key. Managed identity and service principal are not supported.
:::

**Hybrid workloads**

It is possible to use multiple executors in the same Nextflow pipeline. This feature enables the deployment of hybrid workloads in which some tasks are executed locally (or on an HPC cluster) and some tasks are offloaded to Azure Batch.

Use process selectors in your Nextflow configuration to execute specific processes with Azure Batch. For example:

```groovy
process {
    withLabel: bigTask {
        executor = 'azurebatch'
        queue = 'my-batch-pool'
        container = 'my/image:tag'
    }
}

bucketDir = 'az://my-container/some/path'

azure {
    storage {
        accountName = '<STORAGE_ACCOUNT_NAME>'
        accountKey = '<STORAGE_ACCOUNT_KEY>'
    }
    batch {
        location = '<LOCATION>'
        accountName = '<BATCH_ACCOUNT_NAME>'
        accountKey = '<BATCH_ACCOUNT_KEY>'
    }
}
```

With the above configuration, processes labeled as `bigTask` run on Azure Batch and use `bucketDir` as the work directory, while all other processes run locally and use the local work directory specified by `workDir`. Nextflow handles file transfers between the local and remote work directories.

Launch the pipeline with the above configuration:

```bash
nextflow run <PIPELINE_NAME>
```

:::{note}
When using Fusion, `bucketDir` is not needed as all tasks will use the remote work directory. Specify `workDir` instead.
:::

## Advanced configuration

See the {ref}`config-azure` configuration scope for the list of available configuration options.
