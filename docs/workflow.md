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

(workflow-process-invocation)=

## Process invocation

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

### Process named outputs

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

### Process named stdout

The `emit` option can also be used to name a `stdout` output:

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

:::{note}
Optional params for a process input/output are always prefixed with a comma, except for `stdout`. Because `stdout` does not have an associated name or value like other types, the first param should not be prefixed.
:::

## Subworkflows

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

### Workflow parameters

A workflow component can access any variable or parameter defined in the global scope:

```groovy
params.data = '/some/data/file'

workflow my_pipeline {
    if( params.data )
        bar(params.data)
    else
        bar(foo())
}
```

### Workflow inputs

A workflow can declare one or more input channels using the `take` keyword. For example:

```groovy
workflow my_pipeline {
    take: data

    main:
    foo(data)
    bar(foo.out)
}
```

Multiple inputs must be specified on separate lines:

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

### Workflow outputs

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

### Workflow named outputs

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

### Workflow entrypoint

A workflow with no name (also known as the *implicit workflow*) is the default entrypoint of the Nextflow pipeline. A different workflow entrypoint can be specified using the `-entry` command line option.

:::{note}
Implicit workflow definitions are ignored when a script is included as a module. This way, a workflow script can be written in such a way that it can be used either as a library module or an application script.
:::

### Workflow composition

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
