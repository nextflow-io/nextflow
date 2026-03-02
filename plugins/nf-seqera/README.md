# Seqera Executor plugin for Nextflow

## Summary

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
    accessToken = '<SEQERA ACCESS TOKEN>'
}

```

Alternatively, set the access token via environment variable:

```bash
export TOWER_ACCESS_TOKEN='<YOUR ACCESS TOKEN>'
```

## Examples

### Running a workflow with the Seqera executor

`nextflow.config`:

```groovy
plugins {
    id 'nf-seqera'
}

process {
    executor = 'seqera'
}

tower {
    accessToken = '<SEQERA ACCESS TOKEN>'
}

seqera {
    executor {
        region = 'eu-west-1'
    }
}
```

`main.nf`:

```groovy
process HELLO {
    output:
    path 'hello.txt'

    script:
    '''
    echo "Hello from Seqera Cloud" > hello.txt
    '''
}

workflow {
    HELLO()
}
```

### Using resource labels for cost tracking

```groovy
seqera {
    executor {
        region = 'us-east-1'
        labels = [team: 'genomics', project: 'wgs-analysis']
        autoLabels = true
    }
}
```

### Using the resource prediction model

```groovy
seqera {
    executor {
        region = 'eu-west-1'
        predictionModel = 'qr/v1'
    }
}
```

## Resources

- [Seqera Platform Documentation](https://docs.seqera.io/)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
