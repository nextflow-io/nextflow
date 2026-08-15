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

#### Processes depending on `task.memory`

When a prediction model is enabled the scheduler can allocate less memory than the task requested.
The task script however is rendered *before* the task is scheduled, therefore a script referencing
`task.memory` carries the memory that was *requested*, not the one that was allocated e.g.

```groovy
process FOO {
    memory 8.GB

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar app.jar
    """
}
```

Here `-Xmx8g` is baked into the command even when the scheduler allocates less, and the task fails
with an out-of-memory error. To prevent this the executor submits the affected tasks with prediction
model `none` and reports a warning.

The references are collected at compile time from the process AST, so the check also covers the
directives that are resolved independently of the task script, and the values defined in the config
file — including the common `ext.args` idiom:

```groovy
process { withName: FOO { ext.args = { "-Xmx${task.memory.toGiga()}g" } } }
```

Set the `seqera/predictionModel` hint explicitly on the process to override this behaviour:

```groovy
process FOO {
    hints 'seqera/predictionModel': 'qr/v1'
}
```

Note that `task.cpus` is not subject to this check, and neither are the `template` files, whose
content is only read at run time.

## Resources

- [Seqera Platform Documentation](https://docs.seqera.io/)

## License

[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)
