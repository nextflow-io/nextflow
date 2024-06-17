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

## Implicit workflow

A script can define a single workflow without a name (also known as the *implicit workflow*), which is the default entrypoint of the script. The `-entry` command line option can be used to execute a different workflow as the entrypoint at runtime.

:::{note}
Implicit workflow definitions are ignored when a script is included as a module. This way, a script can be written such that it can be either imported as a module or executed as a pipeline.
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
The use of global variables and params in named workflows is discouraged because it breaks the modularity of the workflow. As a best practice, every workflow input should be explicitly defined as such in the `take:` block, and params should only be used in the implicit workflow.
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

:::{note}
This feature requires the `nextflow.preview.output` feature flag to be enabled.
:::

A script may define the set of outputs that should be published by the implicit workflow, known as the workflow output definition:

```groovy
workflow {
    foo(bar())
}

output {
    directory 'results'
}
```

The output definition must be defined after the implicit workflow.

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
