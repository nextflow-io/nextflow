(workflow-page)=

# Workflows

In Nextflow, a **workflow** is a function that is specialized for composing processes and dataflow logic (i.e. channels and operators).

See {ref}`syntax-workflow` for a full description of the workflow syntax.

:::{note}
Workflows were introduced in DSL2. If you are still using DSL1, see {ref}`dsl1-page` for more information about how to migrate your Nextflow pipelines to DSL2.
:::

## Entry workflow

A script can define up to one *entry workflow*, which does not have a name and serves as the entrypoint of the script:

```nextflow
workflow {
    channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
        .map { v -> "$v world!" }
        .view()
}
```

### Parameters

Parameters can be declared by assigning a `params` property to a default value:

```nextflow
params.input = '/some/data/file'
params.save_intermeds = false

workflow {
    if( params.input )
        analyze(params.input, params.save_intermeds)
    else
        analyze(fake_input(), params.save_intermeds)
}
```

:::{note}
As a best practice, params should be used only in the entry workflow and passed to workflows and processes as explicit inputs.
:::

The default value can be overridden by the command line, params file, or config file. Parameters from multiple sources are resolved in the order described in {ref}`cli-params`.

## Named workflows

A *named workflow* is a workflow that can be called by other workflows:

```nextflow
workflow my_workflow {
    hello()
    bye( hello.out.collect() )
}

workflow {
    my_workflow()
}
```

The above example defines a workflow named `my_workflow` which is called by the entry workflow. Both `hello` and `bye` could be any other process or workflow.

### Takes and emits

The `take:` section is used to declare the inputs of a named workflow:

```nextflow
workflow my_workflow {
    take:
    data1
    data2

    main:
    hello(data1, data2)
    bye(hello.out)
}
```

Inputs can be specified like arguments when calling the workflow:

```nextflow
workflow {
    my_workflow( channel.of('/some/data') )
}
```

The `emit:` section is used to declare the outputs of a named workflow:

```nextflow
workflow my_workflow {
    main:
    hello(data)
    bye(hello.out)

    emit:
    bye.out
}
```

When calling the workflow, the output can be accessed using the `out` property, i.e. `my_workflow.out`.

If an output is assigned to a name, the name can be used to reference the output from the calling workflow. For example:

```nextflow
workflow my_workflow {
    main:
    hello(data)
    bye(hello.out)

    emit:
    my_data = bye.out
}
```

The result of the above workflow can be accessed using `my_workflow.out.my_data`.

:::{note}
Every output must be assigned to a name when multiple outputs are declared.
:::

(workflow-process-invocation)=

## Calling processes and workflows

Processes and workflows are called like functions, passing their inputs as arguments:

```nextflow
process hello {
    output:
    path 'hello.txt', emit: txt

    script:
    """
    your_command > hello.txt
    """
}

process bye {
    input:
    path 'hello.txt'

    output:
    path 'bye.txt', emit: txt

    script:
    """
    another_command hello.txt > bye.txt
    """
}

workflow hello_bye {
    take:
    data

    main:
    hello()
    bye(data)
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    hello_bye(data)
}
```

Processes and workflows have a few extra rules for how they can be called:

- Processes and workflows can only be called by workflows

- A given process or workflow can only be called once in a given workflow. To use a process or workflow multiple times in the same workflow, use {ref}`module-aliases`.

The "return value" of a process or workflow call is the process outputs or workflow emits, respectively. The return value can be assigned to a variable or passed into another call:

```nextflow
workflow hello_bye {
    take:
    data

    main:
    bye_out = bye(hello(data))

    emit:
    bye_out
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    bye_out = hello_bye(data)
}
```

Named outputs can be accessed as properties of the return value:

```nextflow
workflow hello_bye {
    take:
    data

    main:
    hello_out = hello(data)
    bye_out = bye(hello_out.txt)

    emit:
    bye = bye_out.txt
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    flow_out = hello_bye(data)
    bye_out = flow_out.bye
}
```

As a convenience, process and workflow outputs can also be accessed without first assigning to a variable, by using the `.out` property of the process or workflow name:

```nextflow
workflow hello_bye {
    take:
    data

    main:
    hello(data)
    bye(hello.out)

    emit:
    bye = bye.out
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    hello_bye(data)
    hello_bye.out.bye.view()
}
```

:::{note}
Process named outputs are defined using the `emit` option on a process output. See {ref}`naming process outputs <process-naming-outputs>` for more information.
:::

:::{note}
Process and workflow outputs can also be accessed by index (e.g., `hello.out[0]`, `hello.out[1]`, etc.). As a best practice, multiple outputs should be accessed by name.
:::

Workflows can be composed in the same way:

```nextflow
workflow flow1 {
    take:
    data

    main:
    tick(data)
    tack(tick.out)

    emit:
    tack.out
}

workflow flow2 {
    take:
    data

    main:
    tick(data)
    tock(tick.out)

    emit:
    tock.out
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    flow1(data)
    flow2(flow1.out)
}
```

:::{note}
The same process can be called in different workflows without using an alias, like `tick` in the above example, which is used in both `flow1` and `flow2`. The workflow call stack determines the *fully qualified process name*, which is used to distinguish the different process calls, i.e. `flow1:tick` and `flow2:tick` in the above example.
:::

:::{tip}
The fully qualified process name can be used as a {ref}`process selector <config-process-selectors>` in a Nextflow configuration file, and it takes priority over the simple process name.
:::

## Special operators

The following operators have a special meaning when used in a workflow with process and workflow calls.

### Pipe `|`

The `|` *pipe* operator can be used to chain processes, operators, and workflows:

```nextflow
process greet {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

workflow {
    channel.of('Hello','Hola','Ciao')
        | greet
        | map { v -> v.toUpperCase() }
        | view
}
```

The above snippet defines a process named `greet` and invokes it with the input channel. The result is then piped to the {ref}`operator-map` operator, which converts each string to uppercase, and finally to the {ref}`operator-view` operator which prints it.

The same code can also be written as:

```nextflow
workflow {
    ch1 = channel.of('Hello','Hola','Ciao')
    ch2 = greet( ch1 )
    ch2.map { v -> v.toUpperCase() }.view()
}
```

### And `&`

The `&` *and* operator can be used to call multiple processes in parallel with the same channel(s):

```nextflow
process greet {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

process to_upper {
    input:
    val data

    output:
    val result

    exec:
    result = data.toUpperCase()
}

workflow {
    channel.of('Hello')
        | map { v -> v.reverse() }
        | (greet & to_upper)
        | mix
        | view
}
```

In the above snippet, the initial channel is piped to the {ref}`operator-map` operator, which reverses the string value. Then, the result is passed to the processes `greet` and `to_upper`, which are executed in parallel. Each process outputs a channel, and the two channels are combined using the {ref}`operator-mix` operator. Finally, the result is printed using the {ref}`operator-view` operator.

The same code can also be written as:

```nextflow
workflow {
    ch = channel.of('Hello').map { v -> v.reverse() }
    ch_greet = greet(ch)
    ch_upper = to_upper(ch)
    ch_greet.mix(ch_upper).view()
}
```

(workflow-recursion)=

## Process and workflow recursion

:::{versionadded} 21.11.0-edge
:::

:::{note}
This feature requires the `nextflow.preview.recursion` feature flag to be enabled.
:::

Processes can be invoked recursively using the `recurse` method.

```{literalinclude} snippets/recurse-process.nf
:language: nextflow
```

```{literalinclude} snippets/recurse-process.out
:language: console
```

In the above example, the `count_down` process is first invoked with the value `params.start`. On each subsequent iteration, the process is invoked again using the output from the previous iteration. The recursion continues until the specified condition is satisfied, as defined by the `until` method, which terminates the recursion.

The recursive output can also be limited using the `times` method:

```groovy
count_down
    .recurse(params.start)
    .times(3)
    .view { v -> "${v}..." }
```

Workflows can also be invoked recursively:

```{literalinclude} snippets/recurse-workflow.nf
:language: nextflow
```

```{literalinclude} snippets/recurse-workflow.out
:language: console
```

**Limitations**

- A recursive process or workflow must have matching inputs and outputs, such that the outputs for each iteration can be supplied as the inputs for the next iteration.

- Recursive workflows cannot use *reduction* operators such as `collect`, `reduce`, and `toList`, because these operators cause the recursion to hang indefinitely after the initial iteration.

(workflow-output-def)=

## Workflow outputs

:::{versionadded} 24.04.0
:::

:::{versionchanged} 24.10.0
A second preview version was introduced. See the {ref}`migration notes <workflow-outputs-second-preview>` for details.
:::

:::{versionchanged} 25.04.0
A third preview version was introduced. See the {ref}`migration notes <workflow-outputs-third-preview>` for details.
:::

:::{note}
This feature requires the `nextflow.preview.output` feature flag to be enabled.
:::

A script can define an *output block* which declares the top-level outputs of the workflow. Each output should be assigned in the `publish` section of the entry workflow. Any channel in the workflow can be assigned to an output, including process and subworkflow outputs. This approach is intended to replace the {ref}`publishDir <process-publishdir>` directive.

Here is a basic example:

```nextflow
process fetch {
    // ...

    output:
    path 'sample.txt'

    // ...
}

workflow {
    main:
    fetch(params.input)

    publish:
    samples = fetch.out
}

output {
    samples {
        path '.'
    }
}
```

In the above example, the output of process `fetch` is assigned to the `samples` workflow output. How this output is published to a directory structure is described in the next section.

(workflow-publishing-files)=

### Publishing files

The top-level output directory of a workflow run can be set using the `-output-dir` command-line option or the `outputDir` config option:

```bash
nextflow run main.nf -output-dir 'my-results'
```

```groovy
// nextflow.config
outputDir = 'my-results'
```

The default output directory is `results` in the launch directory.

By default, all output files are published to the output directory. Each output in the output block can define where files are published using the `path` directive. For example:

```nextflow
workflow {
    main:
    ch_step1 = step1()
    ch_step2 = step2(ch_step1)

    publish:
    step1 = ch_step1
    step2 = ch_step2
}

output {
    step1 {
        path 'step1'
    }
    step2 {
        path 'step2'
    }
}
```

The following directory structure will be created:

```
results/
└── step1/
    └── ...
└── step2/
    └── ...
```

All files received by an output will be published into the specified directory. Lists and maps are recursively scanned for nested files. For example:

```nextflow
workflow {
    main:
    ch_samples = channel.of(
        [ [id: 'SAMP1'], [ file('1.txt'), file('2.txt') ] ]
    )

    publish:
    samples = ch_samples // 1.txt and 2.txt will be published
}
```

The `path` directive can also be a closure which defines a custom publish path for each channel value:

```nextflow
workflow {
    main:
    ch_samples = channel.of(
        [id: 'SAMP1', fastq_1: file('1.fastq'), fastq_1: file('2.fastq')]
    )

    publish:
    samples = ch_samples
}

output {
    samples {
        path { sample -> "fastq/${sample.id}/" }
    }
}
```

The above example will publish each channel value to a different subdirectory. In this case, each pair of FASTQ files will be published to a subdirectory based on the sample ID.

The closure can even define a different path for each individual file using the `>>` operator:

```nextflow
output {
    samples {
        path { sample ->
            sample.fastq_1 >> "fastq/${sample.id}/"
            sample.fastq_2 >> "fastq/${sample.id}/"
        }
    }
}
```

Each `>>` specifies a *source file* and *publish target*. The source file should be a file or collection of files, and the publish target should be a directory or file name. If the publish target ends with a slash, it is treated as the directory in which source files are published. Otherwise, it is treated as the target filename of a source file. Only files that are published with the `>>` operator are saved to the output directory.

### Index files

Each output can create an index file of the values that were published. An index file preserves the structure of channel values, including metadata, which is simpler than encoding this information with directories and file names. The index file can be a CSV (`.csv`), JSON (`.json`), or YAML (`.yml`, `.yaml`) file. The channel values should be files, lists, or maps.

For example:

```nextflow
workflow {
    main:
    ch_samples = channel.of(
        [id: 1, name: 'sample 1', fastq_1: '1a.fastq', fastq_2: '1b.fastq'],
        [id: 2, name: 'sample 2', fastq_1: '2a.fastq', fastq_2: '2b.fastq'],
        [id: 3, name: 'sample 3', fastq_1: '3a.fastq', fastq_2: '3b.fastq']
    )

    publish:
    samples = ch_samples
}

output {
    samples {
        path 'fastq'
        index {
            path 'samples.csv'
        }
    }
}
```

The above example will write the following CSV file to `results/samples.csv`:

```
"1","sample 1","results/fastq/1a.fastq","results/fastq/1b.fastq"
"2","sample 2","results/fastq/2a.fastq","results/fastq/2b.fastq"
"3","sample 3","results/fastq/3a.fastq","results/fastq/3b.fastq"
```

You can customize the index file with additional directives, for example:

```nextflow
index {
    path 'samples.csv'
    header true
    sep '|'
}
```

This example will produce the following index file:

```
"id"|"name"|"fastq_1"|"fastq_2"
"1"|"sample 1"|"results/fastq/1a.fastq"|"results/fastq/1b.fastq"
"2"|"sample 2"|"results/fastq/2a.fastq"|"results/fastq/2b.fastq"
"3"|"sample 3"|"results/fastq/3a.fastq"|"results/fastq/3b.fastq"
```

See [Output directives](#output-directives) for the list of available index directives.

### Output directives

The following directives are available for each output in the output block:

`index`
: Create an index file which will contain a record of each published value.

  The following directives are available in an index definition:

  `header`
  : When `true`, the keys of the first record are used as the column names (default: `false`). Can also be a list of column names. Only used for CSV files.

  `path`
  : The name of the index file relative to the base output directory (required). Can be a CSV, JSON, or YAML file.

  `sep`
  : The character used to separate values (default: `','`). Only used for CSV files.

`label`
: Specify a label to be applied to every published file. Can be specified multiple times.

`path`
: Specify the publish path relative to the output directory (default: `'.'`). Can be a path, a closure that defines a custom directory for each published value, or a closure that publishes individual files using the `>>` operator.

Additionally, the following options from the {ref}`workflow <config-workflow>` config scope can be specified as directives:
- `contentType`
- `enabled`
- `ignoreErrors`
- `mode`
- `overwrite`
- `storageClass`
- `tags`

For example:

```nextflow
output {
    samples {
        mode 'copy'
    }
}
```
