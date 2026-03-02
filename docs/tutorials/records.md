(migrating-records)=

# Migrating to records

Nextflow 26.04 introduces the ability to use *records* and *record types* in a Nextflow pipeline. This tutorial demonstrates how to migrate to records using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example.

## Overview

Nextflow pipelines need a way to model composite data, such as a *sample* consisting of a unique identifer, some files, and other metadata.

The standard way to model composite data in Nextflow is with a *tuple*:

```nextflow
sample = tuple([id: '1'], file('1_1.fastq'), file('1_2.fastq'))

channel.of(sample).view { meta, fastq_1, fastq_2 ->
    "id: ${meta.id}, fastq_1: ${fastq_1.name}, fastq_2: ${fastq_2.name}"
}
```

Tuples are easy to use, but have several downsides:

- **Brittle data structures**: Tuples must be repeatedly unpacked and repacked whenever they pass through a channel operator. The tuple can instead be passed as a single argument, in which case tuple elements must be accessed by index. Either way, the user has to remember the names and ordering of the tuple elements.

- **Tight coupling between processes and workflow logic**: When calling a process, the tuple elements must match the process inputs with the same ordering. Getting this wrong can lead to runtime errors that are hard to debug.

- **Limited type checking**: The type checker can validate whether an incoming tuple has the right number of elements and the right element types, but it cannot validate the semantic name of each element, because tuple elements don't have names. For example, the type `Tuple<String,Path,Path>` does not tell you what any of the elements refer to, only their types. As a result, there are many kinds of semantic errors that cannot be caught by the type checker.

Records are a new data type for modelling complex data structures, as an alternative to tuples:

- **Meaningful data structures**: Record fields are accessed by name instead of index, allowing you to think about your data in terms of meaningful names instead of arbitrary indices.

- **Flexible workflow logic**: Records are *duck-typed*, which means that a record can be used as an input as long as it meets the minimum requirements of that input. This makes it much easier to keep related data together in a workflow without requiring additional logic to translate between channels and process inputs/outputs.

- **Robust type checking**: The type checker can validate workflow logic much more thoroughly when using records. For example, it can validate whether an incoming record has all field names and types required by a process input.

Records allow you to model your data more precisely as it flows through your pipeline. It allows your pipeline code to focus on what matters, and not worry about things that don't matter, making it easier to read and understand.

:::{note}
While records are intended as an alternative to tuples, Nextflow will continue to support tuples as a first-class data type. All existing code will continue to work.
:::

## Developer tooling

Records and record types work best with the [Nextflow language server](https://github.com/nextflow-io/language-server) and [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow).

:::{tip}
See {ref}`devenv-page` for instructions on how to setup VS Code and the Nextflow extension.
:::

### Type checking

The language server can check your code for type-related errors, which is especially useful when using records. For example, it can validate that a channel of records has all required fields when it is passed as input to a process.

### Automatic migration

The **Convert script to static types** command described in {ref}`migrating-static-types` automatically converts process inputs and outputs to records when the language version is set to 26.04.

For example, the following process:

```nextflow
process ALIGN {
    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.bai'), emit: bai
}
```

Is converted to:

```nextflow
process ALIGN {
    input:
    sample: Record {
        meta: Map
        reads: Path
        index: Path
    }

    output:
    record(
        meta: sample.meta,
        bam: file('*.bam'),
        bai: file('*.bai')
    )
}
```

The outputs must satisfy the following conditions in order to be converted:

- Every output must be a tuple with two components, with the first component being a `val` output.

- Every output must have an emit name.

If these conditions are met, the tuple components are combined into a single record, with the emit name used as the field name for the second component of each tuple.

## Example: rnaseq-nf

This section demonstrates how to migrate a pipeline to records using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example. The completed migration is available in the [preview-26-04](https://github.com/nextflow-io/rnaseq-nf/tree/preview-26-04) branch.

See {ref}`rnaseq-nf-page` for an introduction to the rnaseq-nf pipeline.

:::{note}
This tutorial assumes you are familiar with {ref}`static types <migrating-static-types>` and {ref}`workflow outputs <migrating-workflow-outputs>`. Refer to those tutorials for an introduction to these concepts.
:::

### Migrating pipeline inputs

The rnaseq-nf pipeline takes a CSV file, also known as a *samplesheet*, containing a row for each input sample. This samplesheet is loaded into a channel of tuples:

```nextflow
read_pairs_ch = channel.of(params.reads)
    .flatMap { csv -> csv.splitCsv() }
    .map { row ->
        tuple(row[0], file(row[1]), file(row[2]))
    }
```

To load the samplesheet as a channel of records, replace the `tuple()` call with `record()` as follows:

```nextflow
read_pairs_ch = channel.of(params.reads)
    .flatMap { csv -> csv.splitCsv() }
    .map { row ->
        record(id: row[0], fastq_1: file(row[1]), fastq_2: file(row[2]))
    }
```

### Migrating workflows

The `RNASEQ` workflow has the following definition after converting it to static types:

```nextflow
workflow RNASEQ {
    take:
    reads_ch        : Channel<Tuple<String,Path,Path>>
    transcriptome   : Path

    main:
    index = INDEX(transcriptome)
    fastqc_ch = FASTQC(reads_ch)
    quant_ch = QUANT(reads_ch, index)
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    samples : Channel<Tuple<String,Path,Path>> = samples_ch
    index   : Value<Path> = index
}
```

The `reads_ch` input needs to be refactored as a channel of records. While we could use `Channel<Record>`, it is better to use a *record type* that specifies the required fields. To do this, we need to examine the processes that use `ch_reads` and construct a record type that encapsulates those uses.

The `reads_ch` input is used by `FASTQC` and `QUANT`, which both have the following input declaration:

```nextflow
    input:
    (id, fastq_1, fastq_2): Tuple<String, Path, Path>
```

Therefore, we can construct a record type that models these requirements. Update the `reads_ch` input as follows:

```nextflow
workflow RNASEQ {
    take:
    read_pairs_ch: Channel<Sample>
    transcriptome: Path

    // ...
}

record Sample {
    id      : String
    fastq_1 : Path
    fastq_2 : Path
}
```

The `Sample` type also matches the updated `reads_pair_ch` channel that is passed as input from the entry workflow.

TODO: workflow emits, `join` operator

### Migrating processes

See {ref}`process-typed-page` for an overview of typed process inputs and outputs.

<h4>FASTQC</h4>

The `FASTQ` process is defined with the following inputs and outputs:

```nextflow
process FASTQC {
    tag id
    conda 'bioconda::fastqc=0.12.1'

    input:
    (id, fastq_1, fastq_2): Tuple<String, Path, Path>

    output:
    tuple(id, file("fastqc_${id}_logs"))

    script:
    """
    fastqc.sh "${id}" "${fastq_1} ${fastq_2}"
    """
}
```

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
process FASTQC {
    tag sample.id
    conda 'bioconda::fastqc=0.12.1'

    input:
    sample: Record {
        id: String
        fastq_1: Path
        fastq_2: Path
    }

    output:
    record(
        id: sample.id,
        fastqc: file("fastqc_${sample.id}_logs")
    )

    script:
    """
    fastqc.sh "${sample.id}" "${sample.fastq_1} ${sample.fastq_2}"
    """
}
```

In the above:

- The tuple input is converted to a record input using the type `Record`. The field types are specified alongside the field names.

- Since the record input cannot be destructured like a tuple, you must define a name for the record itself (e.g., `sample`), and you must update all references to tuple inputs (e.g., replace `id` with `sample.id`).

- The tuple output is converted to a record by using the `record()` function and specifying a name for each record field.

- Whereas tuple elements must be specified in a particular order, record fields can be specified in any order. The records supplied by the calling workflow must have the same field names and types as the process definition.

<h4>QUANT</h4>

The `QUANT` process is defined with the following inputs and outputs:

```nextflow
process QUANT {
    tag id
    conda 'bioconda::salmon=1.10.3'

    input:
    (id, fastq_1, fastq_2): Tuple<String, Path, Path>
    index: Path

    output:
    tuple(id, file("quant_${id}"))

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

To migrate this process, rewrite the inputs and outputs as follows:

```nextflow
process QUANT {
    tag sample.id
    conda 'bioconda::salmon=1.10.3'

    input:
    sample: Record {
        id: String
        fastq_1: Path
        fastq_2: Path
    }
    index: Path

    output:
    record(
        id: sample.id,
        quant: file("quant_${sample.id}")
    )

    script:
    """
    salmon quant \
        --threads ${task.cpus} \
        --libType=U \
        -i ${index} \
        -1 ${sample.fastq_1} \
        -2 ${sample.fastq_2} \
        -o quant_${sample.id}
    """
}
```

<h4>MULTIQC</h4>

The `MULTIQC` process does not need to be updated because it does not use tuples.

<h4>INDEX</h4>

The `INDEX` process does not need to be updated because it does not use tuples.

## Additional resources

See the following links to learn more about records and record types:

- {ref}`script-records`
- {ref}`stdlib-types-record`
- {ref}`process-typed-page`
