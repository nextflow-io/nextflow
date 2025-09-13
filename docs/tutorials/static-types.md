(migrating-static-types)=

# Migrating to static types

Nextflow 25.10 introduces the ability to use *static types* in a Nextflow pipeline. This tutorial demonstrates how to migrate to static types using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example.

## Overview

Static types are a way to specify the types of variables in Nextflow code, both to document the code and enable deeper forms of validation. The Nextflow language server can use type annotations to identify type-related errors during development, without needing to run the code.

While Nextflow inherited type annotation from Groovy, types could only be specified for functions and local variables, and not for Nextflow-specific concepts such as processes, workflows, and pipeline parameters. Additionally, the Groovy type system is significantly larger and more complex than what is required for Nextflow pipelines.

Nextflow 25.10 provides a native way to specify types at every level of a pipeline, from a pipeline parameter to a local variable in a process, using the {ref}`standard types <stdlib-types>` in the Nextflow standard library.

## Developer tooling

Static types are most effective when used with the [Nextflow language server](https://github.com/nextflow-io/language-server), especially in combination with the [VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow).

### Type checking

Static type checking is currently only available through the language server as an experimental feature.

### Automatic migration

The VS Code extension provides a command for automatically migrating Nextflow pipelines to static types. Search for and select **Convert pipeline to static types** in the command palette to migrate the current project.

The language server will attempt to convert every legacy process to a typed process by migrating inputs and outputs to the new syntax.

- In cases where the type of an input or output cannot be inferred (e.g. `val` inputs and outputs), the type will be left unspecified, and the language server will report an error for each case. If a process has an nf-core [meta.yml](https://nf-co.re/docs/guidelines/components/modules#documentation), the language server will use it to infer the types of `val` inputs and outputs.

- File inputs (`file` and `path` qualifiers) are inferred as type `Path` or `Set<Path>` based on (1) the `arity` option, if specified, or (2) the stage name, if specified. Review the converted code to ensure that the correct type is used. Use the `arity` option in the legacy syntax to ensure the most accurate results.

- File outputs (`file` and `path` qualifiers) are translated to `file()` or `files()` based on (1) the `arity` option, if specified, or (2) whether the file name is a glob pattern. Review the converted code to ensure that the correct output function is used. Use the `arity` option in the legacy syntax to ensure the most accurate results.

:::{tip}
The tooling for automatic migration to static types is actively being developed. While this page shows how to perform migrations manually, some or all of these migration steps will become automatic through developer tools such as the language server.
:::

## Example: rnaseq-nf

This section demonstrates how to migrate to static types using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example. To view the completed migration, see the [`preview-25-10`](https://github.com/nextflow-io/rnaseq-nf/tree/preview-25-10) branch of the rnaseq-nf repository.

### Initial version

The [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline performs a basic RNAseq analysis on a collection of FASTQ paired-end reads:

```nextflow
workflow {
    reads_ch = channel.fromFilePairs( params.reads, checkIfExists: true, flat: true ) 

    (samples_ch, index) = RNASEQ( params.transcriptome, reads_ch )

    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    MULTIQC( multiqc_files_ch, params.multiqc )
}

workflow RNASEQ {
    take:
    transcriptome
    reads_ch

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(reads_ch)
    quant_ch = QUANT(index, reads_ch)
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    index = index
    samples = samples_ch
}
```

Each process declares inputs and outputs using the legacy syntax:

```nextflow
process FASTQC {
    // ...

    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("fastqc_${id}_logs"), emit: logs

    // ...
}
```

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
"outdir": {
  "type": "string",
  "format": "directory-path",
  "description": "The output directory where the results will be saved",
  "default": "results"
},
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
"multiqc": {
  "type": "string",
  "format": "directory-path",
  "description": "Directory containing the configuration for MultiQC",
  "default": "${projectDir}/multiqc"
}
```

Using this schema, replace the legacy parameters with the equivalent `params` block:

```nextflow
params {
    // The output directory where the results will be saved
    outdir: Path = 'results'

    // The input read-pair files
    reads: String = "${projectDir}/data/ggal/ggal_gut_{1,2}.fq"

    // The input transcriptome file
    transcriptome: Path = "${projectDir}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

    // Directory containing the configuration for MultiQC
    multiqc: Path = "${projectDir}/multiqc"
}
```

See {ref}`workflow-params-def` for more information about the `params` block.

### Migrating workflows

To migrate a named workflow to static types, specify the type of each input in the `take:` section.

```nextflow
workflow RNASEQ {
    take:
    transcriptome   : Path
    reads_ch        : Channel<Tuple<String,Path,Path>>

    // ...
}
```

These types can be inferred by examining how the `RNASEQ` workflow is called in the entry workflow:

```nextflow
workflow {
    reads_ch = channel.fromFilePairs( params.reads, checkIfExists: true, flat: true ) 

    (samples_ch, index) = RNASEQ( params.transcriptome, reads_ch )

    // ...
}
```

- The parameter `params.transcriptome` has type `Path` from the `params` block.

- The channel `reads_ch` has type `Channel<E>`, where `E` is the type of each value in the channel. The `fromFilePairs()` factory with `flat: true` emits tuples with three elements corresponding to the id and file pair, which has type `Tuple<String,Path,Path>`.

:::{note}
Tuple types can be quite long. A future version of Nextflow will introduce record types as an alternative to tuples with a better developer experience, such as more concise type annotations. The features introduced in Nextflow 25.10 are focused on enabling static types for existing code with minimal changes to workflow logic.
:::

You can also specify types for workflow emits:

```nextflow
workflow RNASEQ {
    take:
    transcriptome   : Path
    reads_ch        : Channel<Tuple<String,Path,Path>>

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(reads_ch)
    quant_ch = QUANT(index, reads_ch)
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    index   : Value<Path> = index
    samples : Channel<Tuple<String,Path,Path>> = samples_ch
}
```

- The variable `index` comes from the output of `INDEX`. `INDEX` is called with a regular value (rather than a channel), so it emits a dataflow value of type `Value<V>`, where `V` is the type of the inner value. Each `INDEX` task returns a single index file of type `Path`, thus the fully-specified type of `index` is `Value<Path>`.

- The variable `samples_ch` comes from the joined output of `FASTQC` and `QUANT`. Each process is called with the channel `reads_ch`, so they each emit a channel. Each task returns a tuple with an id and output file, both of which have type `Tuple<String,Path>`. The `join` operator combines tuples based on a matching key, which in this case produces tuples of type `Tuple<String,Path,Path>`.

Types are not required for emits since they can be inferred automatically, but they are still useful as documentation and as a sanity check -- if the declared emit type doesn't match the type of the assigned value, the language server can report an error.

### Migrating processes

See {ref}`process-typed-page` for an overview of typed process inputs and outputs.

:::{note}
While this section shows how to migrate processes manually, you can use the language server as described above to migrate code automatically.
:::

<h4>FASTQC</h4>

The `FASTQ` process is defined with the following inputs and outputs:

```nextflow
process FASTQC {
    // ...

    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("fastqc_${id}_logs"), emit: logs

    // ...
}
```

Migrate the process by rewriting the inputs and outputs as follows:

```nextflow
process FASTQC {
    // ...

    input:
    (id, fastq_1, fastq_2): Tuple<String,Path,Path>

    output:
    logs = tuple(id, file("fastqc_${id}_logs"))

    // ...
}
```

- Inputs of type `Path` are treated like `path` inputs in the legacy syntax. Additionally, if `fastq_1` or `fastq_2` were marked as `Path?`, they could be null, which was not allowed with `path` inputs.

- Outputs are defined as assignments, similar to workflow emits. In this case, `logs` could be omitted since there is only one output.

- Values in the `output:` section can use standard library functions as well as several specialized functions for {ref}`process outputs <process-reference-typed>`. In this case, `tuple()` is the {ref}`standard library function <stdlib-namespaces-global>` (not the `tuple` output qualifier) and `file()` is the process output function (not the standard library function).

:::{note}
The other process sections, such as the directives and the `script:` block, are not shown here because they do not need to be changed. As long as the inputs and outputs declare and reference the same variable names and file patterns, the other process sections will behave the same as before.
:::

<h4>QUANT</h4>

The `QUANT` process is defined with the following inputs and outputs:

```nextflow
process QUANT {
    // ...

    input:
    path index
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("quant_${id}")

    // ...
}
```

Migrate the process by rewriting the inputs and outputs as follows:

```nextflow
process QUANT {
    // ...

    input:
    index   : Path
    (id, fastq_1, fastq_2): Tuple<String,Path,Path>

    output:
    tuple(id, file("quant_${id}"))

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

Migrate the process by rewriting the inputs and outputs as follows:

```nextflow
process MULTIQC {
    // ...

    input:
    logs    : Bag<Path>
    config  : Path

    stage:
    stageAs logs, '*'

    output:
    file('multiqc_report.html')

    // ...
}
```

Since the first `path` input was declared with a file pattern, it requires an explicit *stage directive* to stage the file input under a specific alias. You must also declare a variable name for the input, which in the above example is `logs`. The `stageAs` directive specifies that the value of `logs` should be staged using the glob pattern `*`.

In this case, the stage directive can actually be omitted because staging a file input as `*` is equivalent to the default behavior. Inputs that are {ref}`collections <stdlib-types-iterable>` of files (e.g., `Bag<Path>`) are also staged by default.

:::{note}
In the legacy syntax, the `arity` option can be used to specify whether a `path` qualifier expects a single file or collection of files. When using typed inputs and outputs, this behavior is determined by the type, i.e. `Path` vs `Bag<Path>`.
:::

<h4>INDEX</h4>

The `INDEX` process can be migrated using principles already described in the other processes, so it is left as an exercise for the reader.
