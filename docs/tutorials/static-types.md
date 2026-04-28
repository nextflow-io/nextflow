(migrating-static-types)=

# Migrating to static typing

Nextflow 26.04 brings full support for *static typing* in Nextflow code. This tutorial demonstrates how to migrate to static typing using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example.

:::{note}
Static typing is optional. All existing code will continue to work.
:::

## Overview

Static typing allows you to precisely model and validate the structure of your data as it flows through your pipeline. It consists of several new language features:

- **Type annotations** can be added to inputs and outputs at every level of a pipeline, from pipeline parameters to process inputs and outputs, using the {ref}`standard Nextflow types <stdlib-types>`. These annotations make your code easier to understand and are used by the Nextflow language server to identify type-related errors during development.

- **Records** are a new data structure for modeling composite data. They serve as an alternative to tuples -- whereas tuple elements must be accessed by index, record fields are accessed by name. This allows you to model data with meaningful names instead of keeping track of how tuple elements are ordered.

- **Record types** are custom type definitions that can be used to guarantee a minimum set of requirements for a record in a particular context. Records are *duck-typed*, which means that a record can be used as an input as long as it meets the minimum requirements of that input (given by a record type).

## Developer tooling

Static typing works best with the [Nextflow language server](https://github.com/nextflow-io/language-server) and [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow).

:::{tip}
See {ref}`devenv-page` for instructions on how to setup VS Code and the Nextflow extension.
:::

### Type checking

When using static typing, the language server can check your code for type-related errors. For example, it can validate that a channel of records has all the required fields when it is passed as input to a process.

The language server performs type checking on every script that enables the `nextflow.enable.types` feature flag.

### Automatic migration

The Nextflow VS Code extension provides a command for automatically migrating Nextflow pipelines to static types. To migrate a script, open the Command Palette, search for **Convert script to static types**, and select it.

:::{note}
Automatic migration is an experimental feature and may not be able to convert an entire pipeline to static types. Always review generated code for correctness.
:::

## Example: rnaseq-nf

This section demonstrates how to migrate a pipeline to static typing using [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) as an example. See {ref}`rnaseq-nf-page` for an introduction to the pipeline.

The approach is as follows:

1. Convert legacy parameters to a `params` block
2. Convert the primary input (`params.reads`) from a glob pattern to a samplesheet
3. Convert each process to static typing
4. Convert each workflow to static typing

The completed migration is available in the [preview-26-04](https://github.com/nextflow-io/rnaseq-nf/tree/preview-26-04) branch.

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

See {ref}`workflow-typed-params` for more information about the `params` block.

:::{note}
Parameters used only in the config file should be declared in the config, not in the script. Since rnaseq-nf has no such parameters, all parameters are declared in the script. See {ref}`config-params` for more information.
:::

:::{tip}
The rnaseq-nf pipeline initializes the `reads` and `transcriptome` parameters to a test dataset by default, as it is designed as a toy example. In practice, defaults for test data should be defined in a config profile (e.g., `test`).
:::

(static-types-samplesheet)=

### Loading a samplesheet input

The rnaseq-nf pipeline takes a glob pattern of FASTQ pairs (e.g., `data/ggal/ggal_gut_{1,2}.fq`) and uses the `channel.fromFilePairs()` factory to load the files as a channel of tuples:

```nextflow
read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)
```

Each tuple has three elements -- the sample ID (inferred from the file names) and the two FASTQ files.

This approach will not work with static typing because `fromFilePairs()` does not have a well-defined return type. A more robust way to model a collection of samples is with a *samplesheet*, such as a CSV file specifying samples as rows and sample fields as columns.

Create the following samplesheet to represent the test data:

**`data/allreads.csv`**
```
id,fastq_1,fastq_2
gut,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_gut_1.fq,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_gut_2.fq
liver,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_liver_1.fq,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_liver_2.fq
lung,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_lung_1.fq,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_lung_2.fq
spleen,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_spleen_1.fq,https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/refs/heads/master/data/ggal/ggal_spleen_2.fq
```

Refactor `params.reads` to refer to the samplesheet file path instead of a glob pattern:

```nextflow
params {
    // The input samplesheet of paired-end reads
    reads: Path = "${projectDir}/data/allreads.csv"

    // ...
}
```

Refactor the `read_pairs_ch` to load the samplesheet as a channel of records:

```nextflow
read_pairs_ch = channel.of(params.reads)
    .flatMap { csv -> csv.splitCsv() }
    .map { row ->
        record(id: row[0], fastq_1: file(row[1]), fastq_2: file(row[2]))
    }
```

You can simplify the code further by modeling `params.reads` as a collection of records instead of a file path.

Add a header row to the samplesheet:

```
id,fastq_1,fastq_2
gut,...
liver,...
lung,...
spleen,...
```

Refactor `params.reads` as a collection of records:

```nextflow
params {
    // The input samplesheet of paired-end reads
    reads: List<Sample> = "${projectDir}/data/allreads.csv"

    // ...
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}
```

In the above, `Sample` is a *record type* based on the samplesheet structure. When a file path is supplied to a collection-type parameter (e.g., `List<Sample>`), the file path is automatically loaded and parsed into a collection.

Refactor the `read_pairs_ch` to load the collection into a channel:

```nextflow
read_pairs_ch = channel.fromList(params.reads)
```

:::{note}
Collection-type params can also be loaded from JSON and YAML samplesheets. See {ref}`workflow-typed-params` for more information.
:::

### Migrating processes

See {ref}`process-typed-page` for an overview of typed processes.

:::{note}
You must enable the `nextflow.enable.types` feature flag in each script that uses typed processes.
:::

<h4>FASTQC</h4>

The `FASTQC` process is defined as follows:

```nextflow
process FASTQC {
    tag id
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    path "fastqc_${id}_logs"

    script:
    """
    fastqc.sh "${id}" "${fastq_1} ${fastq_2}"
    """
}
```

To migrate the `FASTQC` process, rewrite the inputs and outputs as follows:

```nextflow
nextflow.enable.types = true

process FASTQC {
    tag id
    conda 'bioconda::fastqc=0.12.1'

    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path
    )

    output:
    record(
        id: id,
        fastqc: file("fastqc_${id}_logs")
    )

    script:
    """
    fastqc.sh "${id}" "${fastq_1} ${fastq_2}"
    """
}
```

The tuple input is converted to a record input using the `record()` destructor. The field types are specified alongside the field names. The `path` input qualifier is replaced by the `Path` type.

Whereas tuple elements must be specified in a particular order, record fields can be specified in any order. The records supplied by the calling workflow must have the same field names and types as the process definition.

The tuple output is converted to a record using the `record()` function and specifying a name for each record field. The `path` output qualifier is replaced by the `file()` function (or `files()` if multiple files are expected). See {ref}`process outputs <process-reference-typed>` for the list of special functions that can be used in the `output:` section to retrieve task outputs.

<h4>QUANT</h4>

The `QUANT` process is defined as follows:

```nextflow
process QUANT {
    tag id
    conda 'bioconda::salmon=1.10.3'

    input:
    tuple val(id), path(fastq_1), path(fastq_2)
    path index

    output:
    path "quant_${id}"

    script:
    """
    salmon quant \
        --threads ${task.cpus} \
        --libType=U \
        -i ${index} \
        -1 ${fastq_1} \
        -2 ${fastq_2} \
        -o quant_${id}
    """
}
```

To migrate the `QUANT` process, rewrite the inputs and outputs as follows:

```nextflow
nextflow.enable.types = true

process QUANT {
    tag id
    conda 'bioconda::salmon=1.10.3'

    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path
    )
    index: Path

    output:
    record(
        id: id,
        quant: file("quant_${id}")
    )

    script:
    """
    salmon quant \
        --threads ${task.cpus} \
        --libType=U \
        -i ${index} \
        -1 ${fastq_1} \
        -2 ${fastq_2} \
        -o quant_${id}
    """
}
```

<h4>MULTIQC</h4>

The `MULTIQC` process is defined as follows:

```nextflow
process MULTIQC {
    conda 'bioconda::multiqc=1.27.1'

    input:
    path '*'
    path config

    output:
    path 'multiqc_report.html'

    script:
    """
    cp ${config}/* .
    echo "custom_logo: \$PWD/nextflow_logo.png" >> multiqc_config.yaml
    multiqc -n multiqc_report.html .
    """
}
```

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
nextflow.enable.types = true

process MULTIQC {
    // ...

    input:
    logs: Set<Path>
    config: Path

    // stage:
    // stageAs logs, '*'

    output:
    file('multiqc_report.html')

    // ...
}
```

In a typed process, file patterns for `path` inputs must be declared using a *stage directive*. In this example, the first input uses the variable name `logs`, and the `stageAs` directive stages the input using the glob pattern `*`.

In this case, you can omit the stage directive because `*` matches Nextflow's default staging behavior. Inputs of type `Path` or a `Path` collection (e.g., `Set<Path>`) are staged by default using the pattern `'*'`.

:::{note}
In a legacy process, you can use the `arity` option to specify whether a `path` qualifier expects a single file or collection of files. When using typed inputs and outputs, the type determines this behavior, i.e., `Path` vs `Set<Path>`.
:::

:::{note}
While `List<Path>` and `Bag<Path>` are also valid path collection types, `Set<Path>` is preferred in this case because it represents an unordered collection of files. You should only use `List<Path>` when you want the collection to be ordered.
:::

<h4>INDEX</h4>

Apply the same migration principles from the previous processes to migrate `INDEX`.

### Migrating workflows

Once you migrate every process called by a workflow to static typing, you can migrate the workflow itself.

See {ref}`workflow-typed-page` for an overview of typed workflows.

:::{note}
You must enable the `nextflow.enable.types` feature flag in each script that uses typed workflows.
:::

<h4>RNASEQ</h4>

The `RNASEQ` workflow is defined as follows:

```nextflow
workflow RNASEQ {
    take:
    read_pairs_ch
    transcriptome

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANT(index, read_pairs_ch)

    emit:
    fastqc = fastqc_ch
    quant = quant_ch
}
```

You can infer the type of each workflow input by examining how the workflow is called. In this case, `RNASEQ` is called by the entry workflow with the following arguments:

```nextflow
workflow {
    read_pairs_ch = channel.fromList(params.reads)

    RNASEQ(read_pairs_ch, params.transcriptome)

    // ...
}
```

You can determine the type of each input as follows:

- The channel `read_pairs_ch` has type `Channel<E>`, where `E` is the type of each value in the channel. It is loaded from `params.reads` which has type `List<Sample>`. Therefore `read_pairs_ch` has type `Channel<Sample>`.

- The parameter `params.transcriptome` has type `Path` as defined in the `params` block.

Specify the workflow input types as follows:

```nextflow
nextflow.enable.types = true

workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Sample>
    transcriptome: Path

    // ...
}
```

The `read_pairs_ch` channel also needs to provide all of the record fields required by downstream processes. It is used by `FASTQC` and `QUANT`, which both declare the following record input:

```nextflow
    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path
    )
```

The `Sample` record type contains all of the required fields.

:::{note}
In this case, the records in `read_pairs_ch` are identical to the record inputs of `FASTQC` and `QUANT`. However, `read_pairs_ch` would still be compatible if it contained additional record fields, as long as it contains the fields required by the two processes.
:::

The `FASTQC` and `QUANT` processes produce the channels `fastqc_ch` and `quant_ch`, both of which have type `Channel<Record>`:

- `fastqc_ch` contains records with the fields `id` and `fastqc`
- `quant_ch` contains records with the fields `id` and `quant`

You can infer this type information from the respective process outputs, as shown in the previous section.

These channels are emitted as the outputs of `RNASEQ`. However, with records it is usually simpler to join related channels into a single channel (e.g., to publish the channel as a {ref}`workflow output <migrating-workflow-outputs>`).

Use the `join` operator to join `fastqc_ch` and `quant_ch` by sample ID:

```nextflow
nextflow.enable.types = true

workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Sample>
    transcriptome: Path

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANT(read_pairs_ch, index)
    samples_ch = fastqc_ch.join(quant_ch, by: 'id')

    // ...
}
```

Finally, the workflow needs to be updated to only emit the `samples_ch` channel. Type annotations are not required for emits, but they are still useful as documentation and as a sanity chcek -- if the declared output type doesn't match the assigned value's type, the language server will report it.

While `samples_ch` could be emitted as type `Channel<Record>`, the best practice to use an explicit record type so that downstream workflows know which record fields are available.

Define a new record type based on the available fields in `samples_ch`:

```nextflow
record AlignedSample {
    id: String
    fastqc: Path
    quant: Path
}
```

Update the workflow to emit `samples_ch` with the new record type:

```nextflow
nextflow.enable.types = true

workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Sample>
    transcriptome: Path

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(read_pairs_ch)
    quant_ch = QUANT(read_pairs_ch, index)
    samples_ch = fastqc_ch.join(quant_ch, by: 'id')

    emit:
    samples: Channel<AlignedSample> = samples_ch
}
```

<h4>Entry workflow</h4>

The entry workflow is defined as follows:

```nextflow
workflow {
    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)

    (fastqc_ch, quant_ch) = RNASEQ(read_pairs_ch, params.transcriptome)

    multiqc_files_ch = fastqc_ch.mix(quant_ch).collect()

    MULTIQC(multiqc_files_ch, params.multiqc)
}
```

Rewrite this workflow based on the updated params, processes, and subworkflows:

```nextflow
nextflow.enable.types = true

workflow {
    read_pairs_ch = channel.fromList(params.reads)

    samples_ch = RNASEQ(read_pairs_ch, params.transcriptome)

    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    MULTIQC(multiqc_files_ch, params.multiqc)
}
```

The `reads` param was refactored as a collection of records, so it is loaded into a channel using `channel.fromList`. It is compatible with the records expected by `RNASEQ`.

The `RNASEQ` workflow now returns a single combined channel, so the `mix` operation is no longer needed. The `flatMap` operator is used to extract the files from each record in `samples_ch`.

(preparing-static-types)=

## Preparing for static typing

While static typing can be adopted progressively with existing code, many coding patterns are not compatible with static typing. Following best practices and avoiding anti-patterns beforehand will make it easier to adopt static typing.

### Use the strict syntax

The {ref}`strict syntax <strict-syntax-page>` is required to use static typing. It is enabled by default in Nextflow 26.04.

Before you migrate to static typing, ensure your code adheres to the strict syntax using `nextflow lint` or the language server.

### Avoid deprecated patterns

When preparing for the strict syntax, try to address {ref}`deprecation warnings <strict-syntax-deprecated>` as much as possible. For example:

```nextflow
Channel.from(1, 2, 3).map { it * 2 }        // deprecated
channel.of(1, 2, 3).map { it -> it * 2 }    // best practice
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

However, `set` and `tap` are not supported in typed workflows. Use standard assignments instead.

### Avoid `|` and `&` dataflow operators

The {ref}`special operators <workflow-special-operators>` `|` and `&` provide shorthands for writing dataflow logic:

```nextflow
channel.of('Hello', 'Hola', 'Ciao')
    | greet
    | map { v -> v.toUpperCase() }
    | view
```

However, these special operators are not supported in typed workflows. Use standard assignments and method calls instead:

```nextflow
ch_input = channel.of('Hello', 'Hola', 'Ciao')
ch_greet = greet(ch_input)
ch_greet
    .map { v -> v.toUpperCase() }
    .view()
```

### Avoid `.out` for process and workflow outputs

The `.out` property can be used to access process and workflow outputs in legacy workflows:

```nextflow
MY_WORKFLOW()
MY_WORKFLOW.out.foo.view()
MY_WORKFLOW.out.bar.view()
```

However, this pattern is not supported in typed workflows. Use standard assignments instead:

```nextflow
my_out = MY_WORKFLOW()
my_out.foo.view()
my_out.bar.view()
```

### Avoid `each` input qualifier

The {ref}`each <process-input-each>` input qualifier is not supported in typed processes. Use the {ref}`operator-combine` operator to create a single tuple channel instead.

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

### Avoid legacy operators

Many {ref}`operators <operator-page>` are not statically typed. While you can still use them in typed workflows, the type checker will not be able to fully validate your code. These operators can usually be replaced by another operator and/or a standard library function.

For example, the `splitCsv` operator is not statically typed. Use `flatMap` and the equivalent {ref}`stdlib-types-path` method instead:

```nextflow
// before
channel.fromPath('samplesheet.csv')
    .splitCsv(sep: ',')
    .view()

// after
channel.fromPath('samplesheet.csv')
    .flatMap { csv -> csv.splitCsv(sep: ',') }
    .view()
```

See {ref}`migrating-static-types-operators` for more information.

## Additional resources

See the following links to learn more about static typing:

- {ref}`process-typed-page`
- {ref}`workflow-typed-page`
- {ref}`stdlib-types`
- {ref}`script-records`
- {ref}`syntax-record-type`
