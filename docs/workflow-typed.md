(workflow-typed-page)=

# Workflows (typed)

Typed workflows use a new syntax for inputs and outputs that supports static typing.

To use this feature:

1. Enable the {ref}`strict syntax <strict-syntax-page>` by setting the `NXF_SYNTAX_PARSER` environment variable to `v2`:

    ```bash
    export NXF_SYNTAX_PARSER=v2
    ```

2. Set `nextflow.preview.types = true` in every script that uses typed workflows. The `params` block and `output` block can be used without this feature flag.

See {ref}`syntax-workflow-typed` for the complete syntax reference and {ref}`migrating-static-types` to migrate existing code to static typing.

(workflow-typed-params)=

## Typed parameters

:::{versionadded} 25.10.0
:::

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

Typed parameters should only be referenced in the entry workflow or `output` block. Parameters can be passed to workflows and processes as explicit inputs.

The default value can be overridden by the command line, params file, or config file. Parameters from multiple sources are resolved in the order described in {ref}`cli-params`. Parameters specified on the command line are converted to the appropriate type based on the corresponding type annotation.

A parameter that doesn't specify a default value is a *required* parameter. If a required parameter is not given a value at runtime, the run will fail.

:::{versionadded} 26.04.0
:::

Boolean parameters that don't specify a default value will default to `false`.

Parameters with a collection type (i.e., `List`, `Set`, or `Bag`) can be supplied a file path instead of a literal collection. The file must be CSV, JSON, or YAML. Nextflow will parse the file contents and assign the resulting collection to the parameter. An error is thrown if the file contents do not match the parameter type.

:::{note}
When supplying a CSV file to a collection parameter, the CSV file must contain a header row and must use a comma (`,`) as the column separator.
:::

## Typed outputs

:::{versionadded} 25.10.0
:::

Workflow outputs can use type annotations:

```nextflow
nextflow.preview.types = true

params {
    input: String
}

process fastqc {
    input:
    reads: Path

    output:
    file('fastqc')

    // ...
}

process summary {
    input:
    logs: Path

    output:
    file('report.html')

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

## Typed workflows

:::{versionadded} 26.04.0
:::

:::{warning}
Typed workflows are a preview feature. The syntax and behavior may change in future releases.
:::

Typed workflows can use type annotations in the `take:` and `emit:` sections:

```nextflow
nextflow.preview.types = true

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

### Restricted syntax

The following syntax patterns are no longer supported in typed workflows:

- Using `Channel` to access channel factories (use `channel` instead)

- Using implicit closure parameters (declare parameters explicitly instead)

- Using `set` or `tap` to assign channels (use standard assignments instead)

- Composing dataflow logic with `|` and `&` (use standard method calls instead)

- Accessing process and workflow outputs via `.out` (use standard assignments instead)

See {ref}`preparing-static-types` for more information.

### Operators

The operator library has been updated to provide first-class support for static typing and records. All operators can be used in both typed workflows and legacy workflows. However, only a {ref}`core subset <operator-typed-page>` of operators are recommended for use with static typing.

See {ref}`migrating-static-types-operators` for best practice guidelines when migrating existing code.
