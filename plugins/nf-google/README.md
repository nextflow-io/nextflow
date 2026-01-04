# Google Cloud plugin for Nextflow

The Google Cloud plugin provides support for Google Cloud Platform (GCP), including Google Cloud Batch as a compute executor and Google Cloud Storage as a file system.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-google'
}
```

Configure your Google Cloud credentials and project:

```groovy
google {
    project = '<YOUR PROJECT ID>'
    location = 'us-central1'
}

process.executor = 'google-batch'
workDir = 'gs://<YOUR BUCKET>/work'
```

Authentication can be done via:
- Application Default Credentials
- Service account JSON key file
- Workload Identity (for GKE)

## Examples

### Basic Google Cloud Batch Configuration

```groovy
plugins {
    id 'nf-google'
}

google {
    project = 'my-gcp-project'
    location = 'europe-west1'
    batch {
        spot = true
    }
}

process.executor = 'google-batch'
workDir = 'gs://my-bucket/work'
```

### Using Service Account

```groovy
google {
    project = 'my-gcp-project'
    location = 'us-central1'
    credentials = '/path/to/service-account.json'
}
```

### Machine Type Configuration

```groovy
process {
    executor = 'google-batch'
    machineType = 'n2-standard-4'
    disk = '100 GB'
}
```

## Resources

- [Google Cloud Batch Executor Documentation](https://nextflow.io/docs/latest/google.html)
- [Google Cloud Storage Documentation](https://nextflow.io/docs/latest/google.html#google-cloud-storage)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
