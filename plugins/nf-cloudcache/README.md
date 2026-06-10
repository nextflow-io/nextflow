# Cloud cache plugin for Nextflow

## Summary

The Cloud cache plugin provides cloud-based caching support for Nextflow pipelines. It enables workflow resume capability when using cloud storage as the work directory.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-cloudcache'
}
```

The plugin is automatically activated when using cloud storage (S3, GS, Azure Blob) as the work directory with resume enabled.

```groovy
workDir = 's3://my-bucket/work'
```

Run your pipeline with the `-resume` flag:

```bash
nextflow run main.nf -resume
```

## Examples

### AWS S3 Cache

```groovy
plugins {
    id 'nf-amazon'
    id 'nf-cloudcache'
}

workDir = 's3://my-bucket/work'

aws {
    region = 'us-east-1'
}
```

### Google Cloud Storage Cache

```groovy
plugins {
    id 'nf-google'
    id 'nf-cloudcache'
}

workDir = 'gs://my-bucket/work'

google {
    project = 'my-project'
    location = 'us-central1'
}
```

### Azure Blob Storage Cache

```groovy
plugins {
    id 'nf-azure'
    id 'nf-cloudcache'
}

workDir = 'az://my-container/work'

azure {
    storage {
        accountName = 'mystorageaccount'
        accountKey = System.getenv('AZURE_STORAGE_KEY')
    }
}
```

## Resources

- [Nextflow Cache and Resume](https://nextflow.io/docs/latest/cache-and-resume.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
