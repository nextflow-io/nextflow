# Workflow outputs

- Authors: Ben Sherman
- Status: accepted
- Date: 2025-10-20
- Tags: lang, workflows

## Summary

Introduce a unified, dataflow-centric way to declare the top-level outputs of a workflow.

## Problem Statement

In Nextflow DSL1, each process used `publishDir` to copy output files from the work directory to an external location. Nextflow DSL2 inherited this approach, which became a problem as pipelines grew larger and more modular:

- **Mismatch with reusable modules**: Publishing rules often depend on how a process is used in a given pipeline. Setting `publishDir` inside a module process makes the module less reusable, because the publish path and mode are fixed in the process definition. Using process selectors in configuration is verbose and fragile.

- **Fragmented outputs**: Publishing logic is scattered across many module files. There is no single place to see what a pipeline produces or to reason about the output structure.

- **Redundant configuration**: Common settings like the base output directory and publish mode must be repeated in every `publishDir` declaration.

- **Mismatch with channels**: Channels carry both files and structured metadata (e.g., sample IDs, quality flags). The `publishDir` directive matches files with glob patterns and cannot capture metadata unless the metadata happens to be written to a file. As a result, structured, self-describing outputs are hard to produce.

## Goals

- Declare all pipeline outputs in a single location alongside the entry workflow.

- Assign outputs from channels rather than from individual process definitions, decoupling pipeline-specific publishing rules from reusable modules.

- Support dynamic and fine-grained file publishing to match common publishing patterns (e.g. directory per sample, directory per pipeline step).

- Support structured index files (CSV, JSON, YAML) that preserve output files with associated metadata.

- Define publishing behavior (mode, overwrite, storage class, etc.) globally in the config.

- Support type annotations on output declarations for documentation and compile-time validation.

## Non-goals

- Removing support for `publishDir` immediately. It should continue to work without modification, although it may eventually be phased out as users migrate away from it.

- Publishing outputs from processes or named workflows. Only the entry workflow has a `publish:` section.

- Defining a JSON schema for workflow outputs. Schema/spec generation will be explored in the future.

## Decision

Introduce the `output` block for declaring workflow outputs. Each output defines how files are published to the output directory, and the format of the index file (if defined).

Introduce the output directory as a first-class concept in Nextflow, as well as the `workflow.output` config scope for controlling publishing behavior.

## Core Capabilities

### Output definition

Workflow outputs consist of an `output` block, which declares each output, and a `publish:` section in the entry workflow, which assigns a dataflow source (channel or value) to each output:

```groovy
workflow {
    main:
    ch_fastqc = FASTQC(ch_reads)
    ch_report = MULTIQC(ch_fastqc.collect())

    publish:
    fastqc = ch_fastqc
    report = ch_report
}

output {
    fastqc: Channel<Path> {
        path 'fastqc'
    }
    report: Path {
        path '.'
    }
}
```

Every output assigned in `publish:` must be declared in the `output` block, and vice versa. A mismatch is a compile-time error.

Each output declaration can specify a type annotation for documentation and type checking support. Type annotations are optional and do not change runtime behavior. The type checker uses them to validate the `publish:` section and the `path` directive.

### Output directory

The top-level output directory defaults to `results` in the launch directory. It can be overridden from the command line or config file:

```bash
nextflow run main.nf -output-dir my-results
```

```groovy
// nextflow.config
outputDir = 'my-results'
```

All publish paths declared in the `output` block are relative to this directory. Absolute paths are not allowed.

### Static and dynamic publish paths

The `path` directive accepts a string for a fixed path, or a closure for per-value paths:

```groovy
output {
    // static: all files go to results/fastq/
    reads {
        path 'fastq'
    }

    // dynamic: results are organized by sample id
    samples {
        path { sample -> "${sample.id}" }
    }
}
```

Nextflow recursively scans channel values for files, including files nested inside lists, maps, records, and tuples. Files that did not originate from the work directory are not published.

### Fine-grained file publishing with `>>`

Within a `path` closure, individual files can be published to different locations using the `>>` operator. Only files explicitly captured with `>>` are published; other files in the value are ignored.

```groovy
output {
    samples {
        path { sample ->
            sample.fastqc       >> "fastqc/"
            sample.bam          >> (params.save_bams ? "align/" : null)
            sample.bam_index    >> (params.save_bams ? "align/" : null)
        }
    }
}
```

The *publish source* (left-hand side) should be a file or collection of files. The *publish target* (right-hand side) should be a relative path. If the target has a trailing slash, then the source is published *into* the target directory; otherwise the source is published *as* the target name.

A `null` target suppresses publishing for that file, and a `null` source is also a no-op. Either way, individual files can be excluded by setting the record field to `null` in workflow logic or by using a param in the publish statement.

### Index files

Each output can generate a structured index file that records each published channel value along with its metadata. Supported formats are CSV, JSON, and YAML.

```groovy
output {
    samples {
        path 'fastq'
        index {
            path 'samples.csv'
            header true
        }
    }
}
```

The index file is a *samplesheet*. It preserves the structure of files and metadata in the published channel, and can be passed as input to downstream pipelines. Metadata fields (sample IDs, quality flags, etc.) do not need to be written to a separate metadata file or encoded into file paths.

Files that did not originate from the work directory are not published, but are still included in the index.

### Global defaults via configuration

Common publish settings can be set globally under the `workflow.output` config scope:

```groovy
// nextflow.config
workflow {
    output {
        mode = 'copy'
        overwrite = 'lenient'
    }
}
```

These defaults can be overridden per-output in the `output` block:

```groovy
// main.nf
output {
    fastqc {
        mode = 'symlink'
        overwrite = true
    }
}
```

## Alternatives

### Publishing from processes and subworkflows

Earlier iterations allowed for workflow outputs to be published from subworkflows or processes, instead of requiring all workflow outputs to be propagated up to the entry workflow.

While this approach is less verbose, it breaks the modularity of processes and subworkflows. Publishing behavior belongs to the pipeline, not to individual subcomponents that could be shared across many pipelines. The process or subworkflow should expose all of its outputs as channels, and the calling pipeline should decide whether and how to publish these outputs.

On the other hand, propagating all workflow outputs to the top will make pipelines more verbose, especially when using "skinny tuple" channels. Migrating from tuples to records reduces this verbosity, so large pipelines should be migrated to records before they are migrated to workflow outputs.

## Links

- Community issues: [#4042](https://github.com/nextflow-io/nextflow/issues/4042), [#4661](https://github.com/nextflow-io/nextflow/issues/4661), [#4670](https://github.com/nextflow-io/nextflow/issues/4670)
- [Workflow params ADR](./20250825-workflow-params.md)
- [Record types ADR](./20260306-record-types.md)
