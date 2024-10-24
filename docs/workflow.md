(workflow-page)=

# Workflows

In Nextflow, a **workflow** is a function that is specialized for composing processes and dataflow logic (i.e. channels and operators).

A script can define up to one *entry workflow*, which does not have a name and serves as the entrypoint of the script:

```groovy
workflow {
    Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
        | map { v -> "$v world!" }
        | view
}
```

A *named workflow*, on the other hand, is a workflow that can be called from other workflows:

```groovy
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

```groovy
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

```groovy
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

```groovy
workflow {
    my_workflow( Channel.of('/some/data') )
}
```

## Workflow outputs (`emit`)

The `emit:` section is used to declare workflow outputs:

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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

```groovy
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
    Channel.of('Hello','Hola','Ciao')
        | foo
        | map { v -> v.toUpperCase() }
        | view
}
```

The above snippet defines a process named `foo` and invokes it with the input channel. The result is then piped to the {ref}`operator-map` operator, which converts each string to uppercase, and finally to the {ref}`operator-view` operator which prints it.

The same code can also be written as:

```groovy
workflow {
    ch1 = Channel.of('Hello','Hola','Ciao')
    ch2 = foo( ch1 )
    ch2.map { v -> v.toUpperCase() }.view()
}
```

### And `&`

The `&` *and* operator can be used to call multiple processes in parallel with the same channel(s):

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
    Channel.of('Hello')
        | map { v -> v.reverse() }
        | (foo & bar)
        | mix
        | view
}
```

In the above snippet, the initial channel is piped to the {ref}`operator-map` operator, which reverses the string value. Then, the result is passed to the processes `foo` and `bar`, which are executed in parallel. Each process outputs a channel, and the two channels are combined using the {ref}`operator-mix` operator. Finally, the result is printed using the {ref}`operator-view` operator.

The same code can also be written as:

```groovy
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

:::{note}
This feature requires the `nextflow.preview.output` feature flag to be enabled.
:::

A script may define the set of outputs that should be published by the entry workflow, known as the workflow output definition:

```groovy
workflow {
    foo(bar())
}

output {
    directory 'results'
}
```

The output definition must be defined after the entry workflow.

### Publishing channels

Processes and workflows can each define a `publish` section which maps channels to publish targets. For example:

```groovy
process foo {
    // ...

    output:
    path 'result.txt', emit: results

    publish:
    results >> 'foo'

    // ...
}

workflow foobar {
    main:
    foo(data)
    bar(foo.out)

    publish:
    foo.out >> 'foobar/foo'

    emit:
    bar.out
}
```

In the above example, the output `results` of process `foo` is published to the target `foo/` by default. However, when the workflow `foobar` invokes process `foo`, it publishes `foo.out` (i.e. `foo.out.results`) to the target `foobar/foo/`, overriding the default target defined by `foo`.

In a process, any output with an `emit` name can be published. In a workflow, any channel defined in the workflow, including process and subworkflow outputs, can be published.

:::{note}
If the publish source is a process/workflow output (e.g. `foo.out`) with multiple channels, each channel will be published. Individual output channels can also be published by index or name (e.g. `foo.out[0]` or `foo.out.results`).
:::

As shown in the example, workflows can override the publish targets of process and subworkflow outputs. This way, each process and workflow can define some sensible defaults for publishing, which can be overridden by calling workflows as needed.

By default, all files emitted by the channel will be published into the specified directory. If a channel emits list values, any files in the list (including nested lists) will also be published. For example:

```groovy
workflow {
    ch_samples = Channel.of(
        [ [id: 'sample1'], file('sample1.txt') ]
    )

    publish:
    ch_samples >> 'samples' // sample1.txt will be published
}
```

### Publish directory

The `directory` statement is used to set the top-level publish directory of the workflow:

```groovy
output {
    directory 'results'

    // ...
}
```

It is optional, and it defaults to the launch directory (`workflow.launchDir`). Published files will be saved within this directory.

### Publish targets

A publish target is a name with a specific publish configuration. By default, when a channel is published to a target in the `publish:` section of a process or workflow, the target name is used as the publish path.

For example, given the following output definition:

```groovy
workflow {
    ch_foo = foo()
    ch_bar = bar(ch_foo)

    publish:
    ch_foo >> 'foo'
    ch_bar >> 'bar'
}

output {
    directory 'results'
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

:::{note}
The trailing slash in the target name is not required; it is only used to denote that the target name is intended to be used as the publish path.
:::

:::{warning}
The target name must not begin with a slash (`/`), it should be a relative path name.
:::

Workflows can also disable publishing for specific channels by redirecting them to `null`:

```groovy
workflow {
    ch_foo = foo()

    publish:
    ch_foo >> (params.save_foo ? 'foo' : null)
}
```

Publish targets can be customized in the output definition using a set of options similar to the {ref}`process-publishdir` directive.

For example:

```groovy
output {
    directory 'results'
    mode 'copy'

    'foo' {
        mode 'link'
    }
}
```

In this example, all files will be copied by default, and files published to `foo/` will be hard-linked, overriding the default option.

Available options:

`contentType`
: *Currently only supported for S3.*
: Specify the media type a.k.a. [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_Types) of published files (default: `false`). Can be a string (e.g. `'text/html'`), or `true` to infer the content type from the file extension.

`enabled`
: Enable or disable publishing (default: `true`).

`ignoreErrors`
: When `true`, the workflow will not fail if a file can't be published for some reason (default: `false`).

`mode`
: The file publishing method (default: `'symlink'`). The following options are available:

  `'copy'`
  : Copy each file into the output directory.

  `'copyNoFollow'`
  : Copy each file into the output directory without following symlinks, i.e. only the link is copied.

  `'link'`
  : Create a hard link in the output directory for each file.

  `'move'`
  : Move each file into the output directory.
  : Should only be used for files which are not used by downstream processes in the workflow.

  `'rellink'`
  : Create a relative symbolic link in the output directory for each file.

  `'symlink'`
  : Create an absolute symbolic link in the output directory for each output file.

`overwrite`
: When `true` any existing file in the specified folder will be overwritten (default: `'standard'`). The following options are available:

  `false`
  : Never overwrite existing files.

  `true`
  : Always overwrite existing files.

  `'deep'`
  : Overwrite existing files when the file content is different.

  `'lenient'`
  : Overwrite existing files when the file size is different.

  `'standard'`
  : Overwrite existing files when the file size or last modified timestamp is different.

`path`
: Specify the publish path relative to the output directory (default: the target name). Can only be specified within a target definition.

`storageClass`
: *Currently only supported for S3.*
: Specify the storage class for published files.

`tags`
: *Currently only supported for S3.*
: Specify arbitrary tags for published files. For example:
  ```groovy
  tags FOO: 'hello', BAR: 'world'
  ```

### Index files

A publish target can create an index file of the values that were published. An index file is a useful way to save the metadata associated with files, and is more flexible than encoding metadata in the file path. Currently only CSV files are supported.

For example:

```groovy
workflow {
    ch_foo = Channel.of(
        [id: 1, name: 'foo 1'],
        [id: 2, name: 'foo 2'],
        [id: 3, name: 'foo 3']
    )

    publish:
    ch_foo >> 'foo'
}

output {
    directory 'results'

    'foo' {
        index {
            path 'index.csv'
        }
    }
}
```

The above example will write the following CSV file to `results/foo/index.csv`:

```csv
"id","name"
"1","foo 1"
"2","foo 2"
"3","foo 3"
```

You can customize the index file by specifying options in a block, for example:

```groovy
index {
    path 'index.csv'
    header ['name', 'extra_option']
    sep '\t'
    mapper { val -> val + [extra_option: 'bar'] }
}
```

The following options are available:

`header`
: When `true`, the keys of the first record are used as the column names (default: `false`). Can also be a list of column names.

`mapper`
: Closure which defines how to transform each published value into a CSV record. The closure should return a list or map. By default, no transformation is applied.

`path`
: The name of the index file relative to the target path (required).

`sep`
: The character used to separate values (default: `','`).
