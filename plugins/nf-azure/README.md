# Azure plugin for Nextflow

This plugin implements the support for Azure Blob storage as file system
provider (via JSR203 interface) and Azure Batch executor for Nextflow.

## Development

Build Nextflow as usual:

```bash
make compile
```

Use the following Nextflow configuration:

```groovy
plugins {
  id 'nf-azure'
}

azure {
  storage {
    accountKey = "<YOUR STORAGE ACCOUNT KEY>"
    accountName = "<YOUR STORAGE ACCOUNT KEY>"
  }

  batch {
    endpoint = 'https://<YOUR BATCH ACCOUNT NAME>.westeurope.batch.azure.com'
    accountName = '<YOUR BATCH ACCOUNT NAME>'
    accountKey = '<YOUR BATCH ACCOUNT KEY>'
  }
}

process.executor = 'azurebatch'
workDir = 'az://<YOUR DATA CONTAINER>/work'
```

Then test the local build as usual:

```bash
./launch.sh run -c nextflow.config rnaseq-nf
```

## Todo

* Currently, the Blob storage service uses NettyHttpClient and Batch service
  uses OkHttp client, duplicating the number of required libraries. In principle
  the Blob service can use OkHttp, adding the following deps, however using that
  Nextflow hangs during the shutdown, apparently because the connection pool used
  by the blob service is not closed timely.

  ```groovy
  compile('com.azure:azure-storage-blob:12.9.0') {
      exclude group: 'org.slf4j', module: 'slf4j-api'
      exclude group: 'com.azure', module: 'azure-core-http-netty'
  }
  compile('com.azure:azure-core-http-okhttp:1.3.3') {
      exclude group: 'org.slf4j', module: 'slf4j-api'
  }
  ```

* Remove invalid directory from .command.run PATH for project having `bin/` folder
* Add the configuration for the region
* Make the backend endpoint optional

### Additional Resources

* https://github.com/Azure/azure-sdk-for-java/wiki
* https://github.com/Azure/azure-sdk-for-java/tree/master/sdk/storage/azure-storage-blob-nio
* https://github.com/Azure/azure-sdk-for-java/blob/master/sdk/storage/azure-storage-blob-nio/src/samples/java/com/azure/storage/blob/nio/ReadmeSamples.java
