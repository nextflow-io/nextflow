# Microsoft Azure plugin for Nextflow

This plugin provides support for Azure Blob Storage as a file system and Azure Batch as a compute executor for Nextflow pipelines.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-azure'
}
```

Configure your Azure credentials and services:

```groovy
azure {
    storage {
        accountName = '<YOUR STORAGE ACCOUNT NAME>'
        accountKey = '<YOUR STORAGE ACCOUNT KEY>'
    }

    batch {
        endpoint = 'https://<YOUR BATCH ACCOUNT NAME>.<REGION>.batch.azure.com'
        accountName = '<YOUR BATCH ACCOUNT NAME>'
        accountKey = '<YOUR BATCH ACCOUNT KEY>'
    }
}
```

Set the executor and work directory:

```groovy
process.executor = 'azurebatch'
workDir = 'az://<YOUR CONTAINER>/work'
```

## Examples

### Basic Azure Batch Configuration

```groovy
plugins {
    id 'nf-azure'
}

azure {
    storage {
        accountName = 'mystorageaccount'
        accountKey = System.getenv('AZURE_STORAGE_KEY')
    }

    batch {
        endpoint = 'https://mybatchaccount.westeurope.batch.azure.com'
        accountName = 'mybatchaccount'
        accountKey = System.getenv('AZURE_BATCH_KEY')
        autoPoolMode = true
        deletePoolsOnCompletion = true
    }
}

process.executor = 'azurebatch'
workDir = 'az://mycontainer/work'
```

### Using Managed Identity

```groovy
azure {
    managedIdentity {
        clientId = '<YOUR MANAGED IDENTITY CLIENT ID>'
    }
}
```

## Resources

- [Azure Batch Executor Documentation](https://nextflow.io/docs/latest/azure.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
