(migrating-workflow-outputs)=

# Migrating to workflow outputs

The {ref}`workflow output definition <workflow-output-def>` is a new way to define the top-level outputs of a workflow. It is a replacement for the {ref}`publishDir <process-publishdir>` directive. This tutorial describes these changes and explains how to migrate from `publishDir` to workflow outputs using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example.

## Overview

In Nextflow DSL1, pipelines were defined in a single script and there was no concept of workflows. Each process used the `publishDir` directive to publish task outputs, which captured output files with glob patterns and copied them from the work directory to an external location.

Nextflow DSL2 introduced workflows and modules, making it easier to develop large and complex pipelines. However, DSL2 retained the same process-based publishing syntax and became unwieldy for several reasons:

- **Mismatch with reusable modules**: Publishing rules often depend on how a process is used in a pipeline. This made it impractical to set `publishDir` in a reusable way for processes that are shared across many pipelines. Publishing rules could be defined in the configuration, but this approach requires extensive use of {ref}`process selectors <config-process-selectors>`, which are difficult to use for large pipelines.

- **Fragmented outputs**: It is difficult to get a concise view of a workflow's outputs when publishing rules are separated across many different modules.

- **Redundant configuration**: Certain settings, such as the base output directory and publish mode, must be repeated for each `publishDir` declaration, leading to duplicated code.

- **Mismatch with channels**: Channels, the primary data structure in Nextflow, contain files and associated metadata that can be accessed by name. However, `publishDir` uses glob patterns to match files, and cannot publish metadata unless it happens to be in a file. This mismatch makes it difficult to translate channels into pipeline outputs.

Workflow outputs were introduced to address these problems by providing a unified, structured, and flexible way to publish outputs:

- **Unified output definition**: Workflow outputs are declared in an `output` block alongside the entry workflow, ensuring that there is a single comprehensive view of what a pipeline produces.

- **Channel-based publishing**: Instead of publishing files from individual processes, workflow outputs are assigned from channels in the entry workflow. The channel itself can be saved as an *index file*, such as a CSV or JSON file, which provides a structured view of the output directory and can be ingested by downstream pipelines.

- **Flexible file selection**: By default, all files in a published channel are included. However, the published channel can be configured to publish specific files by name, instead of using glob patterns. This approach to publishing files is a natural extension of workflows and channels.

- **Simple configuration**: The base output directory is defined as a global configuration setting, and all files are published into this directory. Publish settings such as the mode are also defined as configuration settings under the `workflow.output` scope, reducing code duplication.

## Timeline

Nextflow {ref}`24.04 <workflow-outputs-first-preview>` introduced the workflow output definition as a preview feature. It has since undergone multiple revisions, with the third preview currently available in Nextflow {ref}`25.04 <workflow-outputs-third-preview>`.

Nextflow 25.10 will finalize workflow outputs and bring them out of preview. The `publishDir` directive will continue to be supported, but will be deprecated. It may be removed in a future release.

## Example: rnaseq-nf

This section describes how to migrate from `publishDir` to workflow outputs using the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline as an example. To view the completed migration, see the [`preview-25-04`](https://github.com/nextflow-io/rnaseq-nf/tree/preview-25-04) branch of the rnaseq-nf repository.

See {ref}`rnaseq-nf-page` for an introduction to the rnaseq-nf pipeline.

### Replacing `publishDir` with workflow outputs

Start by removing each `publishDir` directive and publishing the corresponding process output channel in the entry workflow.

Declare an output for each channel in the `output` block and publish the corresponding channel in the `publish:` section of the entry workflow:

```nextflow
workflow {
    main:
    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)

    (fastqc_ch, quant_ch) = RNASEQ(read_pairs_ch, params.transcriptome)

    multiqc_files_ch = fastqc_ch.mix(quant_ch).collect()

    multiqc_report = MULTIQC(multiqc_files_ch, params.multiqc)

    publish:
    fastqc_logs = fastqc_ch
    multiqc_report = multiqc_report
}

output {
    fastqc_logs {
    }

    multiqc_report {
    }
}
```

:::{note}
Each output assigned in the `publish:` section must be declared in the `output` block, and vice versa.
:::

Nextflow copies all files in the published channels into the output directory, which is `results` by default. You can set the output directory using the `outputDir` config setting or the `-output-dir` command-line option.

You can set the publish mode in the config. For example:

```groovy
workflow.output.mode = 'copy'
```

Run the pipeline with the `all-reads` profile to verify the published outputs:

```console
$ nextflow run . -profile conda,all-reads
```

### Customizing the publish paths

The pipeline runs `FASTQC` and `QUANT` for each input sample. However, the workflow publishes only the `FASTQC` results. The workflow passes the `QUANT` results to `MULTIQC` but doesn't publish them directly.

Improve the workflow outputs by also publishing the outputs of `QUANT`:

```nextflow
workflow {
    main:
    // ...

    publish:
    fastqc_logs = fastqc_ch
    quant = quant_ch
    multiqc_report = multiqc_report
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

This directory will quickly become cluttered as you process more samples. It would be better to group the `FASTQC` and `QUANT` results into separate subdirectories:

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

Achieve this directory structure by customizing the `output` block.

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

Configure the `fastqc_logs` and `quant` outputs in the `output` block to use dynamic publish paths:

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

Preserving the sample ID in each output channel allows you to customize the publish path without trying to parse the file name. The dynamic path is applied to each channel value to determine the target name for the given file.

:::{note}
The closure parameters for the dynamic publish path must match the structure of the published channel.
:::

### Generating an index file

An *index file* is a manifest, or *index*, of the published files and their metadata for a workflow output. Nextflow can create an index file for each workflow output by saving the channel as a CSV, JSON, or YAML file.

For example, if you enable the index file for `fastqc_logs`:

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

The workflow produces the following index file:

```console
$ cat results/fastqc.csv
"id","fastqc"
"lung","results/fastqc/lung"
"spleen","results/fastqc/spleen"
"gut","results/fastqc/gut"
"liver","results/fastqc/liver"
```

The index file mirrors the structure of the published channel, and it provides a structured view of the output directory. Index files are equivalent to samplesheets, and can be used as inputs to downstream pipelines.

You could define two index files for `fastqc_logs` and `quant`. However, since these outputs essentially provide different *slices* of data for the same set of samples, you can also combine them into a single output with one index file.

Use the `join` operator to combine the `FASTQC` and `QUANT` results into a single channel:

```nextflow
workflow {
    main:
    // ...

    samples_ch = fastqc_ch
        .join(quant_ch)
        .map { id, fastqc, quant ->
            [id: id, fastqc: fastqc, quant: quant]
        }

    multiqc_files_ch = samples_ch
        .flatMap { sample -> [sample.fastqc, sample.quant] }
        .collect()
    multiqc_report = MULTIQC( multiqc_files_ch, params.multiqc )

    publish:
    samples = samples_ch
    multiqc_report = multiqc_report
}
```

This example uses maps instead of tuples so that you can access fields by name, and so that the index file can use the map keys as column names.

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

Since each channel value now contains multiple files that go to different subdirectories, you must use *publish statements* in the `path` directive to route each file to the appropriate location.

Run the pipeline, then verify the index file:

```console
$ nextflow run . -profile conda,all-reads -resume
$ cat results/samples.csv
"id","fastqc","quant"
"lung","results/fastqc/lung","results/quant/lung"
"gut","results/fastqc/gut","results/quant/gut"
"liver","results/fastqc/liver","results/quant/liver"
"spleen","results/fastqc/spleen","results/quant/spleen"
```

In the future, if you add a tool with per-sample outputs, you only need to join the tool output into the `samples_ch` channel and update the output `path` directive accordingly. This approach keeps the output definition concise as you add more tools to the pipeline. Additionally, a single unified index file for all per-sample outputs is easier for downstream pipelines to consume, rather than cross-referencing multiple related index files.
