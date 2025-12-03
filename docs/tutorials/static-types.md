(migrating-static-types)=

# Migrating to static types

Nextflow 25.10 introduces the ability to use *static types* in a Nextflow pipeline. This tutorial demonstrates how to migrate to static types using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example.

:::{note}
Static types in Nextflow 25.10 are optional. All existing code will continue to work.
:::

## Overview

Static types allow you to specify the types of variables and parameters in Nextflow code. This language feature serves two purposes:

1. **Better documentation**: Type annotations serve as code documentation, making your code more precise and easier to understand.

2. **Better validation**: The Nextflow language server uses these type annotations to perform *type checking*, which allows it to identify type-related errors during development without requiring code execution.

While Nextflow inherited type annotations from Groovy, types were limited to functions and local variables and weren't supported by Nextflow-specific constructs, such as processes, workflows, and pipeline parameters. Additionally, the Groovy type system is significantly larger and more complex than necessary for Nextflow pipelines.

Nextflow 25.10 provides a native way to specify types at every level of a pipeline, from pipeline parameters to process inputs and outputs, using the {ref}`standard types <stdlib-types>` in the Nextflow standard library.

## Developer tooling

Static types work best with the [Nextflow language server](https://github.com/nextflow-io/language-server) and [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow).

:::{tip}
See {ref}`devenv-page` for instructions on how to setup VS Code and the Nextflow extension.
:::

### Type checking

Static type checking is currently available through the language server as an experimental feature.

### Automatic migration

The Nextflow VS Code extension provides a command for automatically migrating Nextflow pipelines to static types. To migrate a script, open the Command Palette, search for **Convert script to static types**, and select it.

:::{note}
The extension can also convert an entire project using the **Convert pipeline to static types** command.
:::

The conversion consists of the following steps:

- Legacy parameter declarations are converted to a `params` block. If a `nextflow_schema.json` file is present, it is used to infer parameter types.

- Legacy processes are converted to typed processes.

The language server uses the following rules when attempting to convert a legacy process:

- When input or output types cannot be inferred (e.g., `val` inputs and outputs), the type remains unspecified and the language server reports an error. However, if a process includes an nf-core [meta.yml](https://nf-co.re/docs/guidelines/components/modules#documentation), the language server uses it to infer the appropriate type.

- File inputs (`file` and `path` qualifiers) are converted to `Path` or `Set<Path>` based on (1) the `arity` option or (2) stage name when specified. To ensure accurate conversion, specify the `arity` option in the legacy syntax and review the converted code to verify the correct type is used.

- File outputs (`file` and `path` qualifiers) are converted to `file()` or `files()` based on the `arity` option when specified. To ensure accurate conversion, specify the `arity` option in the legacy syntax and review the converted code to verify the correct output function is used.

## Example: rnaseq-nf

This section demonstrates how to migrate a pipeline to static types using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example. The completed migration is available in the [preview-25-10](https://github.com/nextflow-io/rnaseq-nf/tree/preview-25-10) branch.

See {ref}`rnaseq-nf-page` for an introduction to the rnaseq-nf pipeline.

:::{tip}
While much of this migration can be performed automatically by the language server, this example is intended to provide a concrete example of a pipeline that uses static types. Always review generated code for correctness.
:::

### Migrating pipeline parameters

The pipeline defines the following parameters in the main script using the legacy syntax:

```nextflow
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"
```

The pipeline also has a `nextflow_schema.json` schema with the following properties:

```json
"reads": {
  "type": "string",
  "description": "The input read-pair files",
  "default": "${projectDir}/data/ggal/ggal_gut_{1,2}.fq"
},
"transcriptome": {
  "type": "string",
  "format": "file-path",
  "description": "The input transcriptome file",
  "default": "${projectDir}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
},
"outdir": {
  "type": "string",
  "format": "directory-path",
  "description": "The output directory where the results will be saved",
  "default": "results"
},
"multiqc": {
  "type": "string",
  "format": "directory-path",
  "description": "Directory containing the configuration for MultiQC",
  "default": "${projectDir}/multiqc"
}
```

To migrate the pipeline parameters, use the schema and legacy parameters to define the equivalent `params` block:

```nextflow
params {
    // The input read-pair files
    reads: String = "${projectDir}/data/ggal/ggal_gut_{1,2}.fq"

    // The input transcriptome file
    transcriptome: Path = "${projectDir}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

    // The output directory where the results will be saved
    outdir: Path = 'results'

    // Directory containing the configuration for MultiQC
    multiqc: Path = "${projectDir}/multiqc"
}
```

See {ref}`workflow-params-def` for more information about the `params` block.

:::{note}
Parameters that are only used in the config file should be declared in the config file, not the script. Since rnaseq-nf does not have any such parameters, all parameters are declared in the script. See {ref}`config-params` for more information.
:::

:::{tip}
The rnaseq-nf pipeline initializes the `reads` and `transcriptome` parameters to a test dataset by default because it is designed to be run as a toy example. In practice, these defaults should be defined in a config profile (e.g., `test`).
:::

### Migrating workflows

The type of each workflow input can be inferred by examining how the workflow is called.

```nextflow
workflow RNASEQ {
    take:
    read_pairs_ch
    transcriptome

    main:
    // ...

    emit:
    fastqc = fastqc_ch
    quant = quant_ch
}
```

The `RNASEQ` workflow is called by the entry workflow with the following arguments:

```nextflow
workflow {
    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)

    (fastqc_ch, quant_ch) = RNASEQ(read_pairs_ch, params.transcriptome)

    // ...
}
```

You can determine the type of each input as follows:

- The channel `read_pairs_ch` has type `Channel<E>`, where `E` is the type of each value in the channel. The `fromFilePairs()` factory with `flat: true` emits tuples containing a sample ID and two file paths. Therefore the type of `read_pairs_ch` is `Channel<Tuple<String,Path,Path>>`.

- The parameter `params.transcriptome` has type `Path` as defined in the `params` block.

To update the workflow inputs, specify the input types as follows:

```nextflow
workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Tuple<String,Path,Path>>
    transcriptome: Path

    // ...
}
```

:::{note}
Type annotations can become long for tuples with many elements. Future Nextflow versions will introduce alternative data types (i.e., record types) that have more concise annotations and provide a better overall experience with static types. The focus of Nextflow 25.10 is to enable type checking for existing code with minimal changes.
:::

You can also specify types for workflow outputs. This is typically not required for type checking, since the output type can usually be inferred from the input types and workflow logic. However, explicit type annotations are still useful as a sanity check -- if the declared output type doesn't match the assigned value's type, the language server can report an error.

You can use the same approach for inputs to determine the type of each output:

- The variable `fastqc_ch` comes from the output of `FASTQC`. This process is called with a channel (`read_pairs_ch`), so it also emits a channel. Each task returns a tuple containing an id and an output file, so the element type is `Tuple<String,Path>`. Therefore, the type of `fastqc_ch` is `Channel<Tuple<String,Path>>`.

- The variable `quant_ch` comes from the output of `QUANT`. THis process is called with a channel (`read_pairs_ch`), so it also emits a channel. Each task returns a tuple containing an id and an output file, so the element type is `Tuple<String,Path>`. Therefore, the type of `quant_ch` is `Channel<Tuple<String,Path>>`.

To update the workflow outputs, specify the output types as follows:

```nextflow
workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Tuple<String,Path,Path>>
    transcriptome: Path

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANT(index, read_pairs_ch)

    emit:
    fastqc: Channel<Tuple<String,Path>> = fastqc_ch
    quant: Channel<Tuple<String,Path>> = quant_ch
}
```

### Migrating processes

See {ref}`process-typed-page` for an overview of typed process inputs and outputs.

<h4>FASTQC</h4>

The `FASTQ` process is defined with the following inputs and outputs:

```nextflow
process FASTQC {
    // ...

    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    path "fastqc_${id}_logs"

    // ...
}
```

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
process FASTQC {
    // ...

    input:
    (id, fastq_1, fastq_2): Tuple<String,Path,Path>

    output:
    file("fastqc_${id}_logs")

    // ...
}
```

In the above:

- Inputs of type `Path` are treated like `path` inputs in the legacy syntax. Additionally, if `fastq_1` or `fastq_2` were marked as `Path?`, they could be null, which is not supported by `path` inputs.

- Outputs are normally defined as assignments, similar to workflow emits. In this case, however, the name was omitted since there is only one output.

- Values in the `output:` section can use several specialized functions for {ref}`process outputs <process-reference-typed>`. Here, `file()` is the process output function (not the standard library function).

:::{note}
Other process sections, such as the directives and the `script:` block, are not shown because they do not require changes. As long as the inputs and outputs declare and reference the same variable names and file patterns, the other process sections will behave the same as before.
:::

<h4>QUANT</h4>

The `QUANT` process is defined with the following inputs and outputs:

```nextflow
process QUANT {
    // ...

    input:
    tuple val(id), path(fastq_1), path(fastq_2)
    path index

    output:
    path "quant_${id}"

    // ...
}
```

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
process QUANT {
    // ...

    input:
    (id, fastq_1, fastq_2): Tuple<String,Path,Path>
    index: Path

    output:
    file("quant_${id}")

    // ...
}
```

<h4>MULTIQC</h4>

The `MULTIQC` process is defined with the following inputs and outputs:

```nextflow
process MULTIQC {
    // ...

    input:
    path '*'
    path config

    output:
    path 'multiqc_report.html'

    // ...
}
```

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
process MULTIQC {
    // ...

    input:
    logs: Set<Path>
    config: Path

    stage:
    stageAs logs, '*'

    output:
    file('multiqc_report.html')

    // ...
}
```

In a typed process, file patterns for `path` inputs must be declared using a *stage directive*. In this example, the first input uses the variable name `logs`, and the `stageAs` directive stages the input using the glob pattern `*`.

In this case, you can omit the stage directive because `*` matches Nextflow's default staging behavior. Inputs of type `Path` or a `Path` collection (e.g., `Set<Path>`) are staged by default using the pattern `'*'`.

:::{note}
In the legacy syntax, you can use the `arity` option to specify whether a `path` qualifier expects a single file or collection of files. When using typed inputs and outputs, the type determines this behavior, i.e., `Path` vs `Set<Path>`.
:::

:::{note}
While `List<Path>` and `Bag<Path>` are also valid path collection types, `Set<Path>` is recommended because it most accurately represents an unordered collection of files. You should only use `List<Path>` when you want the collection to be ordered.
:::

<h4>INDEX</h4>

Apply the same migration principles from the previous processes to migrate `INDEX`.

(preparing-static-types)=

## Preparing for static types

While static types can be adopted progressively with existing code, some coding patterns cannot be statically typed and will prevent the type checker from validating your entire pipeline.

This section provides best practices for writing Nextflow code that works well with static types.

### Use the strict syntax

The {ref}`strict syntax <strict-syntax-page>` remains optional in Nextflow 25.10. However, it is required to use type annotations.

Before you migrate to static types, ensure your code adheres to the strict syntax using `nextflow lint` or the language server.

### Avoid deprecated patterns

When preparing for the strict syntax, try to address {ref}`deprecation warnings <strict-syntax-deprecated>` as much as possible. For example:

```nextflow
Channel.from(1, 2, 3).map { it * 2 }    // deprecated
channel.of(1, 2, 3).map { v -> v * 2 }  // best practice
```

The above example shows how to avoid three deprecated patterns:

1. Using `Channel` to access channel factories (use `channel` instead)
2. Using the deprecated `channel.from` factory (use `channel.of` or `channel.fromList` instead)
3. Using the implicit `it` closure parameter (declare the parameter explicitly instead) 

### Avoid `set` and `tap` operators

Nextflow provides three ways to assign a channel: a standard assignment, the `set` operator, and the `tap` operator:

```nextflow
ch = channel.of(1, 2, 3)            // standard assignment
channel.of(10, 20, 30).set { ch }   // set
channel.of(10, 20, 30).tap { ch }   // tap
```

However, the type checker cannot validate code that uses `set` or `tap`. Use standard assignments in your dataflow logic to enable full type checking.

### Avoid `|` and `&` operators

The {ref}`special operators <workflow-special-operators>` `|` and `&` provide shorthands for writing dataflow logic:

```nextflow
channel.of('Hello', 'Hola', 'Ciao')
    | greet
    | map { v -> v.toUpperCase() }
    | view
```

However, the type checker cannot validate code that uses these special operators. Use standard assignments and method calls instead:

```nextflow
ch_input = channel.of('Hello', 'Hola', 'Ciao')
ch_greet = greet(ch_input)
ch_greet
    .map { v -> v.toUpperCase() }
    .view()
```

### Avoid `each` input qualifier

The {ref}`each <process-input-each>` input qualifier is not supported by typed processes. Use the {ref}`operator-combine` operator to create a single tuple channel instead.

For example:

```nextflow
process align {
    input:
    path seq
    each mode

    script:
    """
    t_coffee -in $seq -mode $mode > result
    """
}

workflow {
    sequences = channel.fromPath('*.fa')
    methods = ['regular', 'espresso', 'psicoffee']

    align(sequences, methods)
}
```

Rewrite the script to use the `combine` operator. It becomes:

```nextflow
process align {
    input:
    tuple path(seq), val(mode)

    script:
    """
    t_coffee -in $seq -mode $mode > result
    """
}

workflow {
    sequences = channel.fromPath('*.fa')
    methods = ['regular', 'espresso', 'psicoffee']

    align(sequences.combine(methods))
}
```

:::{tip}
The `each` qualifier is discouraged in modern Nextflow code. While it provides a convenient shorthand for combining multiple inputs, it couples the process definition with external workflow logic. Since the introduction of DSL2, Nextflow aims to treat processes as standalone modules that are decoupled from workflow logic.
:::

### Avoid functions with ambiguous return types

Some {ref}`standard library <stdlib-types>` functions and {ref}`operators <operator-page>` have ambiguous return types. While you can still use them with static types, the type checker will not be able to fully validate your code.

You can tell whether a function has an ambiguous return type by inspecting the function signature. The return type is ambiguous if it is `?` or contains `?` as a type argument.

For example, the `flatten` operator has the following signature:

```nextflow
def flatten() -> Channel<?>
```

Use alternatives that have well-defined return types. For example, use `flatMap` instead of `flatten`.

:::{tip}
The commonly-used operators listed under {ref}`dataflow-type-channel` have well-defined return types and are recommended for use with static types.
:::

If you cannot avoid an ambiguously-typed function, you can use a *hard cast* to specify the expected type, provided you know what the type will be at runtime.

For example:

```nextflow
channel.fromPath('samplesheet.csv')
    .splitCsv()
    .view()
```

The `splitCsv` operator has an ambiguous return type (`Channel<?>`).  It emits lists by default but emits maps when `header` is enabled, making the return type dependent on runtime conditions.


To address this issue:

1. Replace the `splitCsv` operation with a `flatMap` operation that calls the `splitCsv` method (defined for `Path`).

2. Add a `map` operation that casts the emitted rows as `List<String>`.

It becomes:

```nextflow
channel.fromPath('samplesheet.csv')
    .flatMap { csv -> csv.splitCsv() }
    .map { row -> row as List<String> }
    .view()
```

The resulting channel has type `Channel<List<String>>`, and the type checker can validate any downstream code that uses this channel.

## Additional resources

See the following links to learn more about static types:

- {ref}`process-typed-page`
- {ref}`stdlib-types`
- {ref}`syntax-process-typed`
