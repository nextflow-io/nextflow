(migrating-workflow-outputs)=

# Migrating to workflow outputs

The {ref}`workflow output definition <workflow-output-def>` is a new way to define the top-level outputs of a workflow. It is a replacement for the {ref}`publishDir <process-publishdir>` directive. This guide describes how to migrate from `publishDir` to workflow outputs.

## Overview

In Nextflow DSL1, a pipeline had to be defined entirely in a single script, and there was no concept of workflows. Each process was responsible for *publishing* task outputs as workflow outputs using the `publishDir` directive, which captured output files with glob patterns and copied them from the work directory to an external location.

Nextflow DSL2 introduced workflows and the ability to define pipeline components in reusable modules, which made it easier to write large and complex pipelines. However, DSL2 did not change the way that outputs were published, and as a result, the process-based publishing approach became unwieldy for a number of reasons:

- Because the publishing rules for a process typically depend on the calling pipeline, it is impractical to set `publishDir` in a generic way for a process that is re-used across many pipelines. Publishing rules can instead be defined in the configuration, but this approach requires heavy use of {ref}`process selectors <config-process-selectors>`, which are difficult to use for large pipelines.

- It is difficult to get a concise view of a workflow's outputs when the publishing rules are tied to processes that are separated across many different modules.

- Settings such as the base output directiry and publish mode must be specified for each `publishDir` setting, which leads to a large amount of duplicated code.

- Data in a Nextflow pipeline is typically structured as a channel of files and associated metadata, where each file and metadata field can be accessed by name. However, `publishDir` uses glob patterns to match files, and cannot publish metadata unless it happens to be in a file. This mismatch makes it difficult to translate channels -- the primary data structure in Nextflow -- into pipeline outputs.

The workflow output definition is a new way to publish outputs that addresses the issues described above:

- Workflow outputs are declared in an `output` block alongside the entry workflow, ensuring that there is a single comprehensive view of what a pipeline produces.

- Workflow outputs are assigned by publishing channels in the `publish:` section of the entry workflow, instead of publishing files in processes, allowing more flexibility in publishing.

- When a channel is published, all of the files that it contains are automatically published. Alternatively, the output can select specific files to publish by name, instead of using glob patterns. This design makes it easy to translate channels into pipeline outputs.

- The base output directory is defined as a global configuration setting, and all files are published into this directory. Publish settings such as the mode are also defined as configuration settings under the `workflow.output` scope, reducing code duplication.

- A workflow output can save the published channel as an *index file*, such as a CSV or JSON file, which serves as a *manifest* of the published files and associated metadata. Index files provide a structured view of pipeline outputs and can be easily ingested by downstream pipelines.

## Timeline

The workflow output definition was introduced in Nextflow {ref}`24.04 <workflow-outputs-first-preview>` as a preview feature. It has since undergone multiple revisions, with the third preview currently available in Nextflow {ref}`25.04 <workflow-outputs-third-preview>`.

Workflow outputs will be finalized and brought out of preview in Nextflow 25.10. The `publishDir` directive will continue to be supported, but will be deprecated. At some point in the future, `publishDir` may be removed.

## Example: rnaseq-nf

We will use the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline to demonstrate how to migrate from `publishDir` to workflow outputs. The final version is available on the [`preview-25.04`](https://github.com/nextflow-io/rnaseq-nf/tree/preview-25-04) branch of rnaseq-nf.

### Initial version

This pipeline takes a transcriptome file and a collection of FASTQ samples, and performs a basic RNAseq analysis:

```nextflow
workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true, flat: true ) 
    RNASEQ( params.transcriptome, read_pairs_ch )
    MULTIQC( RNASEQ.out, params.multiqc )
}

workflow RNASEQ {
    take:
    transcriptome
    read_pairs_ch

    main: 
    INDEX(transcriptome)
    FASTQC(read_pairs_ch)
    QUANT(INDEX.out, read_pairs_ch)

    emit: 
    QUANT.out | concat(FASTQC.out) | collect
}
```

The `FASTQC` and `MULTIQC` processes publish output files using `publishDir`:

```nextflow
params.outdir = 'results'

process FASTQC {
    publishDir params.outdir, mode:'copy'

    // ...

    output:
    path "fastqc_${sample_id}_logs"

    // ...
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    // ...

    output:
    path 'multiqc_report.html'

    // ...
}
```

### Replace `publishDir` with workflow outputs

To migrate to workflow outputs, we need to remove the `publishDir` settings and publish the output channels instead. We'll start by simply publishing each channel as a workflow output.

First, emit the `QUANT` and `FASTQC` outputs separately in the `RNASEQ` workflow:

```nextflow
workflow RNASEQ {
    // ...

    emit:
    fastqc = FASTQC.out
    quant = QUANT.out
}
```

Declare an output for each channel in the `output` block and publish each channel in the entry workflow:

```nextflow
workflow {
    main:
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true, flat: true )
    RNASEQ( params.transcriptome, read_pairs_ch )

    multiqc_files_ch = RNASEQ.out.fastqc
        .concat(RNASEQ.out.quant)
        .collect()
    MULTIQC( multiqc_files_ch, params.multiqc )

    publish:
    fastqc_logs = RNASEQ.out.fastqc
    multiqc_report = MULTIQC.out
}

output {
    fastqc_logs {
    }

    multiqc_report {
    }
}
```

All files in the published channels are copied into the output directory, which is `results` by default. It can be set using the `outputDir` config setting or the `-output-dir` command line option. The publish mode can be set in the config:

```groovy
workflow.output.mode = 'copy'
```

Run the pipeline with the `all-reads` profile to verify the published outputs:

```console
$ nextflow run . -profile conda,all-reads
```

### Configure the publish paths

The pipeline runs `FASTQC` and `QUANT` for each input sample, but currently only the `FASTQC` results are published, while the `QUANT` results are passed to `MULTIQC` but are not published directly.

Let's improve the workflow outputs by also publishing the outputs of `QUANT`:

```nextflow
workflow {
    main:
    // ...

    publish:
    fastqc_logs = RNASEQ.out.fastqc
    quant = RNASEQ.out.quant
    multiqc_report = MULTIQC.out
}

output {
    fastqc_logs {
    }

    quant {
    }

    multiqc_report {
    }
}
```

Running the pipeline with the `all-reads` profile will produce the following output directory:

```console
results
├── fastqc_gut_logs
├── fastqc_liver_logs
├── fastqc_lung_logs
├── fastqc_spleen_logs
├── multiqc_report.html
├── quant_gut
├── quant_liver
├── quant_lung
└── quant_spleen
```

This directory will quickly become cluttered as we produce more samples. It would be better to group the `FASTQC` and `QUANT` results into separate subdirectories:

```console
results
├── fastqc
│   ├── gut
│   ├── liver
│   ├── lung
│   └── spleen
├── multiqc_report.html
└── quant
    ├── gut
    ├── liver
    ├── lung
    └── spleen
```

Fortunately, we can achieve this directory structure by customzing the `output` block.

First, update the `FASTQC` and `QUANT` processes to also emit the sample ID alongside the output files:

```nextflow
process FASTQC {
    // ...

    output:
    tuple val(id), path("fastqc_${id}")

    // ...
}

process QUANT {
    // ...

    output:
    tuple val(id), path("quant_${id}")

    // ...
}
```

Update the `output` block to use dynamic publish paths:

```nextflow
output {
    fastqc_logs {
        path { id, fastqc -> "fastqc/${id}" }
    }

    quant {
        path { id, quant -> "quant/${id}" }
    }

    multiqc_report {
    }
}
```

Preserving the sample ID in each output channel allows us to customize the publish path without trying to parse the file name. The dynamic path is applied to each channel value to determine the target name for the given file.

### Generate an index file

Nextflow can create an *index file* for each workflow output by saving the channel as a CSV, JSON, or YAML file.

For example, if we enable the index file for `fastqc_logs`:

```nextflow
output {
    fastqc_logs {
        path { id, fastqc -> "fastqc/${id}" }
        index {
            path 'fastqc.csv'
            header true
        }
    }
}
```

It will produce the following index file:

```console
$ cat results/fastqc.csv
"id","fastqc"
"lung","results/fastqc/lung"
"spleen","results/fastqc/spleen"
"gut","results/fastqc/gut"
"liver","results/fastqc/liver"
```

An index file is like a *manifest* that lists the published files and their metadata. It provides a structured view of the output directory, and it mirrors the structure of the published channel. Index files are also equivalent to samplesheets, making it easy to use the output of one pipeline as an input to another.

We could define two index files for `fastqc_logs` and `quant`, but this approach will become unwieldy as we add more tools to the pipeline. Since these two outputs essentially provide different *slices* of data for the same set of samples, it would be nice to join them into a single output with one index file. This would also make it easier for downstream pipelines to use our outputs, as they can consult a single index file rather than cross-referencing multiple files.

Use the `join` operator to combine the `FASTQC` and `QUANT` results into a single channel:

```nextflow
workflow {
    main:
    // ...

    samples_ch = RNASEQ.out.fastqc
        .join(RNASEQ.out.quant)
        .map { id, fastqc, quant ->
            [id: id, fastqc: fastqc, quant: quant]
        }

    multiqc_files_ch = samples_ch
        .flatMap { sample -> [sample.fastqc, sample.quant] }
        .collect()
    MULTIQC( multiqc_files_ch, params.multiqc )

    publish:
    samples = samples_ch
    multiqc_report = MULTIQC.out
}
```

We use maps instead of tuples so that the map keys can be used as column names in the index file.

Declare the `samples` output with an index file:

```nextflow
output {
    samples {
        path { sample ->
            sample.fastqc >> "fastqc/${sample.id}"
            sample.quant >> "quant/${sample.id}"
        }
        index {
            path 'samples.csv'
            header true
        }
    }

    multiqc_report {
    }
}
```

Since each channel value now has multiple files, we use *publish statements* in the `path` directive to route each file to the appropriate location.

Finally, run the pipeline to verify the index file:

```console
$ nextflow run . -profile conda,all-reads -resume
$ cat results/samples.csv
"id","fastqc","quant"
"lung","results/fastqc/lung","results/quant/lung"
"gut","results/fastqc/gut","results/quant/gut"
"liver","results/fastqc/liver","results/quant/liver"
"spleen","results/fastqc/spleen","results/quant/spleen"
```
