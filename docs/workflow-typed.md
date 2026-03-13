(workflow-typed-page)=

# Workflows (typed)

:::{versionadded} 25.10.0
:::

:::{note}
Static types require the {ref}`strict syntax <strict-syntax-page>`. Set the `NXF_SYNTAX_PARSER` environment variable to `v2` to enable:

```bash
export NXF_SYNTAX_PARSER=v2
```
:::

:::{note}
Typed workflows require the `nextflow.preview.types` feature flag to be enabled in each script that uses them.
:::

Typed workflows provide first-class support for static typing. They can use type annotations to describe workflow inputs and outputs, and they provide a streamlined set of dataflow operators.

(workflow-typed-params)=

## Typed parameters

A script can declare parameters using the `params` block:

```nextflow
params {
    // Path to input data.
    input: Path

    // Whether to save intermediate files.
    save_intermeds: Boolean
}

workflow {
    analyze(params.input, params.save_intermeds)
}
```

All {ref}`standard types <stdlib-types>` except for the dataflow types (`Channel` and `Value`) can be used for parameters.

Tyepd parameters should only be referenced in the entry workflow or `output` block. Parameters can be passed to workflows and processes as explicit inputs.

The default value can be overridden by the command line, params file, or config file. Parameters from multiple sources are resolved in the order described in {ref}`cli-params`. Parameters specified on the command line are converted to the appropriate type based on the corresponding type annotation.

A parameter that doesn't specify a default value is a *required* parameter. If a required parameter is not given a value at runtime, the run will fail. Boolean parameters that don't specify a default value default to `false`.

:::{versionadded} 26.04.0
:::

Parameters with a collection type (i.e., `List`, `Set`, or `Bag`) can be supplied a file path instead of a literal collection. The file must be CSV, JSON, or YAML. Nextflow will parse the file contents and assign the resuling collection to the parameter. An error is thrown if the file contents do not match the parameter type.

:::{note}
When supplying a CSV file to a collection parameter, the CSV file must contain a header row and must use a comma (`,`) as the column separator.
:::

## Typed outputs

Workflow outputs can use type annotations:

```nextflow
params {
    input: String
}

process fastqc {
    input:
    path reads

    output:
    path 'fastqc'

    // ...
}

process summary {
    input:
    path logs

    output:
    path 'report.html'

    // ...
}

workflow {
    main:
    ch_reads = channel.fromPath(params.input)
    ch_fastqc = fastqc(ch_reads)
    summary_report = summary(ch_fastqc.collect())

    publish:
    fastqc = ch_fastqc
    summary_report = summary_report
}

output {
    fastqc: Channel<Path> {
        path '.'
    }

    summary_report: Path {
        path '.'
    }
}
```

In the above example, the workflow declares two outputs:

- `fastqc`: a channel of FastQC results (`Channel<Path>`)
- `summary_report`: a dataflow value containing a summary report (`Value<Path>` or `Path`)

Type annotations are useful for documenting the structure of each workflow output, and they can be used by the language server to validate the type of each published output.

Outputs that receive a dataflow value can be declared as `Value<V>` or `V` for short. In this example, the `summary_report` is declared with type `Path` as a shorthand for `Value<Path>`.

## Named workflows

The takes and emits in a named workflow can use type annotations:

```nextflow
workflow hello_bye {
    take:
    samples: Channel<Path>

    main:
    ch_hello = hello(samples)
    val_bye = bye(ch_hello.collect())

    emit:
    result: Value<Path> = val_bye
}
```

In the above example, `hello_bye` takes a channel of files (`Channel<Path>`) and emits a dataflow value with a single file (`Value<Path>`). See {ref}`stdlib-types` for the list of available types.

(workflow-typed-dataflow)=

## Dataflow

:::{versionadded} 26.04.0
:::

Typed workflows introduce the following new behaviors for dataflow logic:

- Simplified operator library that supports static types and records
- Support for named arguments when calling processes
- Restrictions on legacy syntax patterns

### Operators

Channels in a typed workflow use a new set of operators.

These operators are a subset of the legacy operators. They provide a core set of dataflow operations, as well as first-class support for static types and records.

The `Channel` type supports the following operators:

- `collect`: collect the channel values into a collection
- `cross`: cross product of two channels
- `filter`: emit only the channel values that satisfy a condition
- `flatMap`: emit multiple values for each channel value with a closure
- `groupBy`: group channel values by a grouping key
- `join`: relational join of two channels based on a matching key
- `map`: transform each channel value with a closure
- `mix`: concatenate two channels
- `reduce`: reduce channel values into a single value with an accumulator
- `subscribe`: perform an action for each channel value
- `unique`: emit unique channel values
- `until`: emit each channel value until a stopping condition is satisfied
- `view`: print each channel value

The `Value` type supports the following operators:

- `flatMap`: transform the value into multiple values and emit them as a channel
- `map`: transform the value with a closure
- `mix`: concatenate the value with a channel
- `subscribe`: perform an action with the value
- `view`: print the value

See {ref}`operator-typed-page` for more information about each operator. See {ref}`migrating-typed-operators` to learn how to migrate existing operators to typed workflows.

### Process named arguments

A common pattern is for a workflow to combine a channel with one or more dataflow values into a single input for a process:

```nextflow
nextflow.preview.types = true

process align {
    input:
    input: Record {
        id: String
        fastq: Path
        index: Path
    }

    // ...
}

workflow align_dedup {
    take:
    ch_samples: Channel<Sample>
    index: Value<Path>

    main:
    ch_align_inputs = ch_samples
        .cross(index)
        .map { sample, index -> sample + record(index: index) }

    align( ch_align_inputs )

    // ...
}

record Sample {
    id: String
    fastq: Path
}
```

This pattern requires a `cross` and `map` operation for each dataflow value that must be added, which quickly becomes verbose.

In a typed workflow, the same behavior can be achieved by calling the process with named arguments:

```nextflow
workflow align_dedup {
    take:
    ch_samples: Channel<Sample>
    index: Value<Path>

    main:
    align(ch_samples, index: index)

    // ...
}
```

The named arguments supplied to `align` are automatically added to each record in `ch_samples`, producing the equivalent of `ch_align_inputs` in the previous example.

Named arguments can be used with a process under the following conditions:

- The process declares a single record input  
- The positional argument (i.e. `ch_samples`) is a channel of records  
- The named arguments are regular values or dataflow values, not channels

### Restricted syntax

The following syntax patterns are not supported in typed workflows.

**Using `Channel` to access channel factories**

Channel factories should be accessed using the `channel` namespace instead of the `Channel` type:

```nextflow
Channel.of(1, 2, 3) // incorrect
channel.of(1, 2, 3) // correct
```

See {ref}`stdlib-namespaces-channel` and {ref}`stdlib-types-channel` for more information.

**Implicit closure parameter**

Closures that expect a single parameter should always declare the parameter rather than relying on the implicit `it` parameter:

```nextflow
ch.map { it * 2 }       // incorrect
ch.map { it -> it * 2 } // correct
ch.map { v -> v * 2 }   // correct
```

**Special operators `|` and `&`**

The {ref}`special operators <workflow-special-operators>` `|` and `&` cannot be used. Use standard method calls instead:

```nextflow
// incorrect
channel.of(1, 2, 2, 3)
    | unique
    | view

// correct
channel.of(1, 2, 2, 3)
    .unique()
    .view()
```

**Accessing process and workflow outputs via `.out`**

The implicit `.out` property for process outputs and workflow emits cannot be used. Use standard assignments instead:

```nextflow
// incorrect
MY_WORKFLOW()
MY_WORKFLOW.out.foo.view()
MY_WORKFLOW.out.bar.view()

// correct
my_out = MY_WORKFLOW()
my_out.foo.view()
my_out.bar.view()
```
