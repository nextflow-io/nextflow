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

A named workflow is a "subworkflow" that can be invoked from other workflows. For example:

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
The use of global variables and params in subworkflows is discouraged because it breaks the modularity of the workflow. As a best practice, every workflow input should be explicitly defined as such in the `take:` block, and params should only be used in the implicit workflow.
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

## Invoking subworkflows

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

(workflow-output-dsl)=

## Publishing outputs

:::{versionadded} 24.04.0
:::

A script may define the set of outputs that should be published by the implicit workflow, known as the workflow output definition or "output block":

```groovy
workflow {
    foo(bar())
}

output {
    directory 'results'

    'foo' {
        select foo.out
    }

    'bar' {
        defaults mode: 'copy', pattern: '*.txt'
        select bar.out
    }
}
```

The output block must be defined after the implicit workflow.

### Output directory

The `directory` statement is used to set the top-level output directory of the workflow:

```groovy
output {
    directory 'results'

    // ...
}
```

It is optional, and it defaults to the launch directory (`workflow.launchDir`).

The output directory can also be defined using the `-output-dir` {ref}`command line option <cli-run>` or the `outputDir` {ref}`config option <config-miscellaneous>`.

### Path definitions

Path definitions are used to definte the directory structure of the published outputs. A path definition is a path name followed by a block which defines the outputs to be published within that path. Like directories, path definitions can be nested.

The path name defines a subdirectory within the output directory, or the parent path if the path definition is nested.

For example, given the following output block:

```groovy
output {
    directory 'results'

    'foo' {
        // ...
    }

    'bar' {
        // ...

        'baz' {
            // ...
        }
    }
}
```

The following directory structure will be created by the workflow:

```
results/
└── foo/
    └── ...
└── bar/
    └── baz/
        └── ...
    └── ...
```

The path name may also contain multiple subdirectories separated by a slash `/`:

```groovy
output {
    'foo/bar/baz' {
        // ...
    }
}
```

It is a shorthand for the following:

```groovy
output {
    'foo' {
        'bar' {
            'baz' {
                // ...
            }
        }
    }
}
```

### Channel selectors

The `select` statement is used to select channels to publish:

```groovy
output {
    'foo' {
        select foo.out
    }
}
```

Any channel defined in the implicit workflow can be referenced in a channel selector, including process and subworkflow outputs.

:::{note}
A process/subworkflow output (e.g. `foo.out`) can only be selected directly if it contains a single output channel. Multi-channel outputs must be selected by index or name, e.g. `foo.out[0]` or `foo.out.samples`.
:::

By default, all files emitted by the channel will be published into the specified directory. If a list value emitted by the channel contains any files, including files within nested lists, they will also be published. For example:

```groovy
workflow {
    ch_samples = Channel.of(
        [ [id: 'sample1'], file('sample1.txt') ]
    )
}

output {
    'samples' {
        // sample1.txt will be published
        select ch_samples
    }
}
```

The publishing behavior can be customized further by using [publish options](#publish-options). See that section for more details.

### Topic selectors

:::{note}
This feature requires the `nextflow.preview.topic` feature flag to be enabled.
:::

The `topic` statement can be used to select a channel topic for publishing:

```groovy
output {
    'samples' {
        topic 'samples'

        // equivalent to:
        select Channel.topic('samples')
    }
}
```

Topic selectors are a useful way to select channels which are deeply nested within subworkflows, without needing to propagate them to the top-level workflow. You can use the {ref}`operator-topic` operator or the `topic` option for {ref}`process outputs <process-additional-options>` to send a channel to a given topic.

Like a channel selector, a topic selector publishes every file that it receives by default, and it can specify [publish options](#publish-options).

### Publish options

The publishing behavior can be configured using the same options available in the {ref}`process-publishdir` directive.

There are several ways to define publish options:

- The `directory` statement

- The `defaults` statement, which defines publish options for a path defintion

- Channel and topic selectors

Publish options are resolved in a cascading manner, in which more specific settings take priority.

Consider the following example:

```groovy
output {
    directory 'results', mode: 'copy'

    'samples' {
        select ch_samples, pattern: '*.txt', mode: 'link'

        'md5' {
            defaults mode: 'link'
            // ...
            topic 'md5', mode: 'copy'
        }
    }
}
```

In this example, the following rules are applied:

- All files will be copied by default

- The channel selector `select ch_samples` will publish via hard link, overriding the output directory default. Additionally, only files matching the pattern `*.txt` will be published.

- All files published to `samples/md5` will be hard-linked by default, overriding the output directory default.

- The topic selector `topic 'md5'` will publish via copy, overriding the default from `samples/md5`.

:::{note}
The only option from `publishDir` that is not allowed is `path`, because the publish path is defined using path definitions.
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
