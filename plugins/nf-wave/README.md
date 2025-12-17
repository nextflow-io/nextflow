# Wave containers plugin for Nextflow

The Wave containers plugin provides integration with the Wave container service for dynamic container building, augmentation, and on-demand software provisioning.

## Get started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-wave'
}
```

Enable Wave in your configuration:

```groovy
wave {
    enabled = true
}
```

Wave automatically builds containers from Conda packages or Dockerfiles defined in your processes.

## Examples

### Basic Wave configuration

```groovy
plugins {
    id 'nf-wave'
}

wave {
    enabled = true
}

process {
    conda = 'samtools=1.17'
}
```

### Using Wave with Conda packages

```groovy
wave {
    enabled = true
}

process ALIGN {
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.17'

    script:
    '''
    bwa mem ref.fa reads.fq | samtools view -bS - > aligned.bam
    '''
}
```

### Wave with Fusion File System

```groovy
wave {
    enabled = true
}

fusion {
    enabled = true
}

process.executor = 'awsbatch'
workDir = 's3://my-bucket/work'
```

### Container Augmentation

```groovy
wave {
    enabled = true
    strategy = 'conda,container'
}

process TOOL {
    container 'ubuntu:22.04'
    conda 'bioconda::samtools=1.17'

    script:
    '''
    samtools --version
    '''
}
```

### Using Private Registries

```groovy
wave {
    enabled = true
}

docker {
    registry = 'quay.io'
}
```

## Resources

- [Wave Documentation](https://docs.seqera.io/wave)
- [Wave Integration Guide](https://nextflow.io/docs/latest/wave.html)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
