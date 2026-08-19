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

> NOTE: When a process references `task.memory` in the script (e.g. `-Xmx${task.memory.toGiga()}g`), resource prediction is disabled for that process. Otherwise, the scheduler could reduce the memory allocation, and the task would try to allocate more memory than is available, and fail with an out-of-memory error.

The same reference can be hidden in a dynamic directive, where it cannot be inspected:

```groovy
process {
    withName: FOO {
        ext.args = { "-Xmx${task.memory.toGiga()}g" }
    }
}
```

Resource prediction is therefore also disabled for a process that defines `ext`, `beforeScript` or `afterScript` with a closure, whether or not the closure actually depends on the task resources. Use a static value, or set the `seqera/predictionModel` hint explicitly, to keep the prediction enabled.

Set the `seqera/predictionModel` hint explicitly on the process to override this behaviour:

```groovy
process {
    withName: FOO {
        hints = ['seqera/predictionModel': 'qr/v1']
    }
}
```

## Resources

- [Seqera Platform Documentation](https://docs.seqera.io/)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
