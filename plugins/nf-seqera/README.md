# Seqera Executor plugin for Nextflow

The Seqera Executor plugin provides integration with Seqera Cloud for executing Nextflow tasks using Seqera's managed compute infrastructure.

## Get Started

To use this plugin, add it to your `nextflow.config`:

```groovy
plugins {
    id 'nf-seqera'
}
```

Configure the Seqera executor:

```groovy
process {
    executor = 'seqera'
}

seqera {
    executor {
        region = 'eu-west-1'
        autoLabels = true
    }
}

tower {
    accessToken = '<YOUR ACCESS TOKEN>'
}

```

Alternatively, set the access token via environment variable:

```bash
export TOWER_ACCESS_TOKEN='<YOUR ACCESS TOKEN>'
```

## Examples

### Basic Configuration

```groovy
plugins {
    id 'nf-seqera'
}

process {
    executor = 'seqera'
}

tower {
    accessToken = '<YOUR ACCESS TOKEN>'
}
```

## Resources

- [Seqera Platform Documentation](https://docs.seqera.io/)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
