(workflow-page)=

# Workflows

In Nextflow, a **workflow** is a function that is specialized for composing processes and dataflow logic (i.e. channels and operators).

A script can define up to one *entry workflow*, which does not have a name and serves as the entrypoint of the script:

```nextflow
workflow {
    Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
        | map { v -> "$v world!" }
        | view
}
```

A *named workflow*, on the other hand, is a workflow that can be called from other workflows:

```nextflow
workflow my_workflow {
    foo()
    bar( foo.out.collect() )
}

workflow {
    my_workflow()
}
```

The above example defines a workflow named `my_workflow` which can be called from another workflow as `my_workflow()`. Both `foo` and `bar` could be any other process or workflow.

See {ref}`syntax-workflow` for a full description of the workflow syntax.

:::{note}
Workflows were introduced in DSL2. If you are still using DSL1, see {ref}`dsl1-page` for more information about how to migrate your Nextflow pipelines to DSL2.
:::

## Using parameters

Parameters can be defined in the script with a default value that can be overridden from the CLI, params file, or config file. Params should only be used by the entry workflow:

```nextflow
params.data = '/some/data/file'

workflow {
    if( params.data )
        bar(params.data)
    else
        bar(foo())
}
```

:::{note}
While params can also be used by named workflows, this practice is discouraged. Named workflows should receive their inputs explicitly through the `take:` section.
:::

## Workflow inputs (`take`)

The `take:` section is used to declare workflow inputs:

```nextflow
workflow my_workflow {
    take:
    data1
    data2

    main:
    foo(data1, data2)
    bar(foo.out)
}
```

Inputs can be specified like arguments when calling the workflow:

```nextflow
workflow {
    my_workflow( Channel.of('/some/data') )
}
```

## Workflow outputs (`emit`)

The `emit:` section is used to declare workflow outputs:

```nextflow
workflow my_workflow {
    main:
    foo(data)
    bar(foo.out)

    emit:
    bar.out
}
```

When calling the workflow, the output can be accessed using the `out` property, i.e. `my_workflow.out`.

If an output is assigned to a name, the name can be used to reference the output from the calling workflow. For example:

```nextflow
workflow my_workflow {
    main:
    foo(data)
    bar(foo.out)

    emit:
    my_data = bar.out
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
process foo {
    output:
    path 'foo.txt', emit: txt

    script:
    """
    your_command > foo.txt
    """
}

process bar {
    input:
    path x

    output:
    path 'bar.txt', emit: txt

    script:
    """
    another_command $x > bar.txt
    """
}

workflow flow {
    take:
    data

    main:
    foo()
    bar(data)
}

workflow {
    data = Channel.fromPath('/some/path/*.txt')
    flow(data)
}
```

Processes and workflows have a few extra rules for how they can be called:

- Processes and workflows can only be called by workflows

- A given process or workflow can only be called once in a given workflow. To use a process or workflow multiple times in the same workflow, use {ref}`module-aliases`.

The "return value" of a process or workflow call is the process outputs or workflow emits, respectively. The return value can be assigned to a variable or passed into another call:

```nextflow
workflow flow {
    take:
    data

    main:
    bar_out = bar(foo(data))

    emit:
    bar_out
}

workflow {
    data = Channel.fromPath('/some/path/*.txt')
    flow_out = flow(data)
}
```

Named outputs can be accessed as properties of the return value:

```nextflow
workflow flow {
    take:
    data

    main:
    foo_out = foo(data)
    bar_out = bar(foo_out.txt)

    emit:
    bar = bar_out.txt
}

workflow {
    data = Channel.fromPath('/some/path/*.txt')
    flow_out = flow(data)
    bar_out = flow_out.bar
}
```

As a convenience, process and workflow outputs can also be accessed without first assigning to a variable, by using the `.out` property of the process or workflow name:

```nextflow
workflow flow {
    take:
    data

    main:
    foo(data)
    bar(foo.out)

    emit:
    bar = bar.out
}

workflow {
    data = Channel.fromPath('/some/path/*.txt')
    flow(data)
    flow.out.bar.view()
}
```

:::{note}
Process named outputs are defined using the `emit` option on a process output. See {ref}`naming process outputs <process-naming-outputs>` for more information.
:::

:::{note}
Process and workflow outputs can also be accessed by index (e.g., `foo.out[0]`, `foo.out[1]`, etc.). Multiple outputs should instead be accessed by name.
:::

Workflows can be composed in the same way:

```nextflow
workflow flow1 {
    take:
    data

    main:
    foo(data)
    bar(foo.out)

    emit:
    bar.out
}

workflow flow2 {
    take:
    data

    main:
    foo(data)
    baz(foo.out)

    emit:
    baz.out
}

workflow {
    data = Channel.fromPath('/some/path/*.txt')
    flow1(data)
    flow2(flow1.out)
}
```

:::{note}
The same process can be called in different workflows without using an alias, like `foo` in the above example, which is used in both `flow1` and `flow2`. The workflow call stack determines the *fully qualified process name*, which is used to distinguish the different process calls, i.e. `flow1:foo` and `flow2:foo` in the above example.
:::

:::{tip}
The fully qualified process name can be used as a {ref}`process selector <config-process-selectors>` in a Nextflow configuration file, and it takes priority over the simple process name.
:::

## Special operators

The following operators have a special meaning when used in a workflow with process and workflow calls.

### Pipe `|`

The `|` *pipe* operator can be used to chain processes, operators, and workflows:

```nextflow
process foo {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

workflow {
    Channel.of('Hello','Hola','Ciao')
        | foo
        | map { v -> v.toUpperCase() }
        | view
}
```

The above snippet defines a process named `foo` and invokes it with the input channel. The result is then piped to the {ref}`operator-map` operator, which converts each string to uppercase, and finally to the {ref}`operator-view` operator which prints it.

The same code can also be written as:

```nextflow
workflow {
    ch1 = Channel.of('Hello','Hola','Ciao')
    ch2 = foo( ch1 )
    ch2.map { v -> v.toUpperCase() }.view()
}
```

### And `&`

The `&` *and* operator can be used to call multiple processes in parallel with the same channel(s):

```nextflow
process foo {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

process bar {
    input:
    val data

    output:
    val result

    exec:
    result = data.toUpperCase()
}

workflow {
    Channel.of('Hello')
        | map { v -> v.reverse() }
        | (foo & bar)
        | mix
        | view
}
```

In the above snippet, the initial channel is piped to the {ref}`operator-map` operator, which reverses the string value. Then, the result is passed to the processes `foo` and `bar`, which are executed in parallel. Each process outputs a channel, and the two channels are combined using the {ref}`operator-mix` operator. Finally, the result is printed using the {ref}`operator-view` operator.

The same code can also be written as:

```nextflow
workflow {
    ch = Channel.of('Hello').map { v -> v.reverse() }
    ch_foo = foo(ch)
    ch_bar = bar(ch)
    ch_foo.mix(ch_bar).view()
}
```

(workflow-output-def)=

## Publishing outputs

:::{versionadded} 24.04.0
:::

:::{versionchanged} 24.10.0
A second preview version has been introduced. Read the [migration notes](#migrating-from-first-preview) for details.
:::

:::{note}
This feature requires the `nextflow.preview.output` feature flag to be enabled.
:::

A workflow can publish outputs by sending channels to "publish targets" in the workflow `publish` section. Any channel in the workflow can be published, including process and subworkflow outputs. This approach is intended to replace the {ref}`publishDir <process-publishdir>` directive.

Here is a basic example:

```nextflow
process foo {
    // ...

    output:
    path 'result.txt', emit: results

    // ...
}

process bar {
    // ...
}

workflow {
    main:
    foo(data)
    bar(foo.out)

    publish:
    foo.out.results >> 'foo'
    bar.out >> 'bar'
}
```

In the above example, the `results` output of process `foo` is published to the target `foo`, and all outputs of process `bar` are published to the target `bar`.

A "publish target" is simply a name that identifies a group of related outputs. How these targets are saved into a directory structure is described in the next section.

:::{tip}
A workflow can override the publish targets of a subworkflow by "re-publishing" the same channels to a different target. However, the best practice is to define all publish targets in the entry workflow, so that all publish targets are defined in one place at the top-level.
:::

### Output directory

The top-level output directory of a workflow run can be set using the `-output-dir` command-line option or the `outputDir` config option:

```bash
nextflow run main.nf -output-dir 'my-results'
```

```groovy
// nextflow.config
outputDir = 'my-results'
```

It defaults to `results` in the launch directory. All published outputs will be saved into this directory.

Each publish target is saved into a subdirectory of the output directory. By default, the target name is used as the directory name.

For example, given the following publish targets:

```nextflow
workflow {
    main:
    ch_foo = foo()
    ch_bar = bar(ch_foo)

    publish:
    ch_foo >> 'foo'
    ch_bar >> 'bar'
}
```

The following directory structure will be created:

```
results/
└── foo/
    └── ...
└── bar/
    └── ...
```

:::{warning}
Target names cannot begin or end with a slash (`/`).
:::

By default, all files emitted by a published channel will be published into the specified directory. If a channel emits list values, each file in the list (including nested lists) will be published. For example:

```nextflow
workflow {
    main:
    ch_samples = Channel.of(
        [ [id: 'foo'], [ file('1.txt'), file('2.txt') ] ]
    )

    publish:
    ch_samples >> 'samples' // 1.txt and 2.txt will be published
}
```

A workflow can also disable publishing for a specific channel by redirecting it to `null`:

```nextflow
workflow {
    main:
    ch_foo = foo()

    publish:
    ch_foo >> (params.save_foo ? 'foo' : null)
}
```

### Customizing outputs

The output directory structure can be customized further in the "output block", which can be defined alongside an entry workflow. The output block consists of "target" blocks, which can be used to customize specific targets.

For example:

```nextflow
workflow {
    // ...
}

output {
    'foo' {
        enabled params.save_foo
        path 'intermediates/foo'
    }

    'bar' {
        mode 'copy'
    }
}
```

This output block has the following effect:

- The target `foo` will be published only if `params.save_foo` is enabled, and it will be published to a different path within the output directory.

- The target `bar` will publish files via copy instead of symlink.

See [Reference](#reference) for all available directives in the output block.

:::{tip}
The output block is only needed if you want to customize the behavior of specific targets. If you are satisfied with the default behavior and don't need to customize anything, the output block can be omitted.
:::

### Dynamic publish path

The `path` directive in a target block can also be a closure which defines a custom publish path for each channel value:

```nextflow
workflow {
    main:
    ch_fastq = Channel.of( [ [id: 'SAMP1'], file('1.fastq'), file('2.fastq') ] )

    publish:
    ch_fastq >> 'fastq'
}

output {
    'fastq' {
        path { meta, fastq_1, fastq_2 -> "fastq/${meta.id}" }
    }
}
```

The above example will publish each channel value to a different subdirectory. In this case, each pair of FASTQ files will be published to a subdirectory based on the sample ID.

The closure can even define a different path for each individual file by returning an inner closure, similar to the `saveAs` option of the {ref}`publishDir <process-publishdir>` directive:

```nextflow
output {
    'fastq' {
        path { meta, fastq_1, fastq_2 ->
            { file -> "fastq/${meta.id}/${file.baseName}" }
        }
    }
}
```

The inner closure will be applied to each file in the channel value, in this case `fastq_1` and `fastq_2`.

:::{tip}
A mapping closure should usually have only one parameter. However, if the incoming values are tuples, the closure can specify a parameter for each tuple element for more convenient access, also known as "destructuring" or "unpacking".
:::

### Index files

A publish target can create an index file of the values that were published. An index file preserves the structure of channel values, including metadata, which is simpler than encoding this information with directories and file names. The index file can be CSV (`.csv`) or JSON (`.json`).

For example:

```nextflow
workflow {
    main:
    ch_fastq = Channel.of(
        [ [id: 1, name: 'sample 1'], '1a.fastq', '1b.fastq' ],
        [ [id: 2, name: 'sample 2'], '2a.fastq', '2b.fastq' ],
        [ [id: 3, name: 'sample 3'], '3a.fastq', '3b.fastq' ]
    )

    publish:
    ch_fastq >> 'fastq'
}

output {
    'fastq' {
        index {
            path 'index.csv'
        }
    }
}
```

The above example will write the following CSV file to `results/fastq/index.csv`:

```csv
"id","name","fastq_1","fastq_2"
"1","sample 1","results/fastq/1a.fastq","results/fastq/1b.fastq"
"2","sample 2","results/fastq/2a.fastq","results/fastq/2b.fastq"
"3","sample 3","results/fastq/3a.fastq","results/fastq/3b.fastq"
```

You can customize the index file with additional directives, for example:

```nextflow
index {
    path 'index.csv'
    header ['id', 'fastq_1', 'fastq_1']
    sep '\t'
    mapper { meta, fq_1, fq_2 -> meta + [fastq_1: fq_1, fastq_2: fq_2] }
}
```

This example will produce the same index file as above, but with the `name` column removed and with tabs instead of commas.

See [Reference](#reference) for the list of available index directives.

### Migrating from first preview

The first preview of workflow publishing was introduced in 24.04. The second preview, introduced in 24.10, made the following breaking changes:

- The process `publish:` section has been removed. Channels should be published only in workflows, ideally the entry workflow.

- The `directory` output directive has been replaced with the `outputDir` config option and `-output-dir` command line option, which is `results` by default. The other directives such as `mode` have been replaced with config options under `workflow.output.*`.

  In other words, only target blocks can be specified in the output block, but target blocks can still specify directives such as `mode`.

- Target names cannot begin or end with a slash (`/`);

### Reference

The following directives are available in a target block:

`index`
: Create an index file which will contain a record of each published value.

  The following directives are available in an index definition:

  `header`
  : When `true`, the keys of the first record are used as the column names (default: `false`). Can also be a list of column names. Only used for `csv` files.

  `mapper`
  : Closure which defines how to transform each published value into a record. The closure should return a list or map. By default, no transformation is applied.

  `path`
  : The name of the index file relative to the target path (required). Can be a `csv` or `json` file.

  `sep`
  : The character used to separate values (default: `','`). Only used for `csv` files.

`path`
: Specify the publish path relative to the output directory (default: the target name). Can be a path, a closure that defines a custom directory for each published value, or a closure that defines a custom path for each individual file.

Additionally, the following options from the {ref}`workflow <config-workflow>` config scope can be specified as directives:
- `contentType`
- `enabled`
- `ignoreErrors`
- `mode`
- `overwrite`
- `storageClass`
- `tags`

:::{note}
Similarly to process directives vs {ref}`process <config-process>` config options, directives in the `output` block are specified without an equals sign (`=`).
:::
