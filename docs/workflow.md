(workflow-page)=

# Workflows

In Nextflow, a **workflow** is a composition of processes and dataflow logic (i.e. channels and operators).

The workflow definition starts with the keyword `workflow`, followed by an optional name, and finally the workflow body delimited by curly braces. A basic workflow looks like the following example:

```groovy
workflow {
    foo()
}
```

Where `foo` could be a function, a process, or another workflow.

Workflows are *lazily executed*, which means that Nextflow parses the entire workflow structure first, and then executes the entire workflow at once. The order in which a task is executed is determined only by its dependencies, so a task will be executed as soon as all of its required inputs are available.

The syntax of a workflow is defined as follows:

```groovy
workflow [ name ] {

    take:
    < workflow inputs >

    main:
    < dataflow statements >

    emit:
    < workflow outputs >

}
```

:::{tip}
The `main:` label can be omitted if there are no `take:` or `emit:` blocks.
:::

:::{note}
Workflows were introduced in DSL2. If you are still using DSL1, see the {ref}`dsl1-page` page to learn how to migrate your Nextflow pipelines to DSL2.
:::

## Entry workflow

A script can define a single workflow without a name (also known as the *entry workflow*), which is the default entrypoint of the script. The `-entry` command line option can be used to execute a different workflow as the entrypoint at runtime.

:::{note}
Entry workflow definitions are ignored when a script is included as a module. This way, a script can be written such that it can be either imported as a module or executed as a pipeline.
:::

## Named workflows

A named workflow is a workflow that can be invoked from other workflows. For example:

```groovy
workflow my_pipeline {
    foo()
    bar( foo.out.collect() )
}

workflow {
    my_pipeline()
}
```

The above snippet defines a workflow named `my_pipeline`, that can be invoked from another workflow as `my_pipeline()`, just like any other function or process.

## Using variables and params

A workflow can access any variable or parameter defined in the global scope:

```groovy
params.data = '/some/data/file'

workflow {
    if( params.data )
        bar(params.data)
    else
        bar(foo())
}
```

:::{tip}
The use of global variables and params in named workflows is discouraged because it breaks the modularity of the workflow. As a best practice, every workflow input should be explicitly defined as such in the `take:` block, and params should only be used in the entry workflow.
:::

## Workflow inputs (`take`)

A workflow can declare one or more input channels using the `take` keyword. For example:

```groovy
workflow my_pipeline {
    take:
    data1
    data2

    main:
    foo(data1, data2)
    bar(foo.out)
}
```

:::{warning}
When the `take` keyword is used, the beginning of the workflow body must be defined with the `main` keyword.
:::

Inputs can be specified like arguments when invoking the workflow:

```groovy
workflow {
    my_pipeline( channel.from('/some/data') )
}
```

## Workflow outputs (`emit`)

A workflow can declare one or more output channels using the `emit` keyword. For example:

```groovy
workflow my_pipeline {
    main:
    foo(data)
    bar(foo.out)

    emit:
    bar.out
}
```

When invoking the workflow, the output channel(s) can be accessed using the `out` property, i.e. `my_pipeline.out`. When multiple output channels are declared, use the array bracket notation or the assignment syntax to access each output channel as described for [process outputs](#process-outputs).

### Named outputs

If an output channel is assigned to an identifier in the `emit` block, the identifier can be used to reference the channel from the calling workflow. For example:

```groovy
workflow my_pipeline {
    main:
    foo(data)
    bar(foo.out)

    emit:
    my_data = bar.out
}
```

The result of the above workflow can be accessed using `my_pipeline.out.my_data`.

(workflow-process-invocation)=

## Invoking processes

A process can be invoked like a function in a workflow definition, passing the expected input channels like function arguments. For example:

```groovy
process foo {
    output:
    path 'foo.txt'

    script:
    """
    your_command > foo.txt
    """
}

process bar {
    input:
    path x

    output:
    path 'bar.txt'

    script:
    """
    another_command $x > bar.txt
    """
}

workflow {
    data = channel.fromPath('/some/path/*.txt')
    foo()
    bar(data)
}
```

:::{warning}
A process can be only be invoked once in a single workflow, however you can get around this restriction by using {ref}`module-aliases`.
:::

### Process composition

Processes with matching input/output declarations can be composed so that the output of the first process is passed as input to the second process. The previous example can be rewritten as follows:

```groovy
workflow {
    bar(foo())
}
```

### Process outputs

A process output can be accessed using the `out` attribute on the corresponding process object. For example:

```groovy
workflow {
    foo()
    bar(foo.out)
    bar.out.view()
}
```

When a process defines multiple output channels, each output can be accessed by index (`out[0]`, `out[1]`, etc.) or by name (see below).

The process output(s) can also be accessed like the return value of a function:

```groovy
workflow {
    f_out = foo()
    (b1, b2) = bar(f_out)
    b1.view()
}
```

#### Named outputs

The `emit` option can be added to the process output definition to assign a name identifier. This name can be used to reference the channel from the calling workflow. For example:

```groovy
process foo {
    output:
    path '*.bam', emit: samples_bam

    '''
    your_command --here
    '''
}

workflow {
    foo()
    foo.out.samples_bam.view()
}
```

When referencing a named output directly from the process invocation, you can use a more concise syntax:

```groovy
workflow {
    ch_samples = foo().samples_bam
}
```

See {ref}`process outputs <process-additional-options>` for more details.

#### Named stdout

The `emit` option can also be used to name a `stdout` output. However, while process output options are usually prefixed with a comma, this is not the case for `stdout`. This is because `stdout` does not have an argument like other types.


```groovy
process sayHello {
    input:
    val cheers

    output:
    stdout emit: verbiage

    script:
    """
    echo -n $cheers
    """
}

workflow {
    things = channel.of('Hello world!', 'Yo, dude!', 'Duck!')
    sayHello(things)
    sayHello.out.verbiage.view()
}
```

## Invoking workflows

Named workflows can be invoked and composed just like any other process or function.

```groovy
workflow flow1 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}

workflow flow2 {
    take: data
    main:
        foo(data)
        baz(foo.out)
    emit:
        baz.out
}

workflow {
    take: data
    main:
        flow1(data)
        flow2(flow1.out)
}
```

:::{note}
Each workflow invocation has its own scope. As a result, the same process can be invoked in two different workflow scopes, like `foo` in the above snippet, which is used in both `flow1` and `flow2`. The workflow execution path, along with the process names, determines the *fully qualified process name* that is used to distinguish the different process invocations, i.e. `flow1:foo` and `flow2:foo` in the above example.
:::

:::{tip}
The fully qualified process name can be used as a {ref}`process selector <config-process-selectors>` in a Nextflow configuration file, and it takes priority over the simple process name.
:::

## Special operators

### Pipe `|`

The `|` *pipe* operator can be used to compose Nextflow processes and operators. For example:

```groovy
process foo {
    input:
    val data

    output:
    val result

    exec:
    result = "$data world"
}

workflow {
   channel.from('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
}
```

The above snippet defines a process named `foo` and invokes it with the `data` channel. The result is then piped to the {ref}`operator-map` operator, which converts each string to uppercase, and finally to the {ref}`operator-view` operator which prints it.

:::{tip}
Statements can also be split across multiple lines for better readability:

```groovy
workflow {
    channel.from('Hello','Hola','Ciao')
      | foo
      | map { it.toUpperCase() }
      | view
}
```
:::

### And `&`

The `&` *and* operator can be used to feed multiple processes with the same channel(s). For example:

```groovy
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
    channel.from('Hello')
      | map { it.reverse() }
      | (foo & bar)
      | mix
      | view
}
```

In the above snippet, the initial channel is piped to the {ref}`operator-map` operator, which reverses the string value. Then, the result is passed to the processes `foo` and `bar`, which are executed in parallel. Each process outputs a channel, and the two channels are combined using the {ref}`operator-mix` operator. Finally, the result is printed using the {ref}`operator-view` operator.

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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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
