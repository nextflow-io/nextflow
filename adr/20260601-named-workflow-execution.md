# Direct execution for named workflows

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-06-01
- Tags: lang, workflows

## Summary

Introduce the ability to execute named workflows directly with the `nextflow run` command.

## Problem Statement

Consider the following entry workflow which simply wraps a named workflow:

```groovy
params {
    input: Path
    index: Path
}

workflow {
    main:
    samples = params.input.splitCsv(header: true) as List<Sample>
    ch_samples = channel.fromList(samples)
    rnaseq = RNASEQ(ch_samples, index)

    publish:
    aligned = rnaseq.aligned
    multiqc_report = rnaseq.multiqc_report
}

output {
    aligned: Channel<AlignedSample> {
        path { s -> /* ... */ }
        index { path 'aligned.json' }
    }
    multiqc_report: Path {
        path '.'
    }
}
```

Where the `RNASEQ` workflow is defined as follows:

```groovy
workflow RNASEQ {
    take:
    ch_samples: Channel<Sample>
    index: Path

    main:
    ch_aligned = ALIGN(ch_samples, index)
    multiqc_report = MULTIQC(ch_aligned.collect())

    emit:
    aligned: Channel<AlignedSample> = ch_aligned
    multiqc_report: Path = multiqc_report
}

record Sample { /* ... */ }
record AlignedSample { /* ... */ }
```

This example demonstrates that most of the `params` / `workflow` / `output` trio can be equivalently expressed by a named workflow: the `params` block mirrors the `take:` section, and the `output` block and `publish:` section together mirror the `emit:` section.

Named workflows typically consume and produce channels so that they can be composed into larger pipelines. But this prevents them from being directly executable -- the purpose of an entry workflow is to translate between dataflow logic and the external world. If this translation could be inferred automatically, it would allow a named workflow to be both executable and composable, eliminating the need to define explicit entry workflows.

## Solution

Given a named workflow with dataflow inputs and outputs, provide a way to execute it directly via `nextflow run`:

- Load each input record channel from an index file (e.g. CSV, JSON, or YAML file)
- Save each output record channel as an index file
- Refer to output files by work directory path instead of publishing them

### Loading input channels from samplesheets

To support direct execution of named workflows, Nextflow needs to load record channels directly from samplesheets. This can be done by loading the samplesheet data based on the file extension (CSV, JSON, YAML), casting each record to the given record type, and loading the collection as a channel.

For example, given the following named workflow:

```groovy
workflow RNASEQ {
    take:
    ch_samples: Channel<Sample>

    // ...
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}
```

Nextflow will generate the following entry workflow:

```groovy
params {
    samples: Path // can be a CSV, JSON, or YAML file
}

workflow {
    ch_samples = channel.fromList(loadData(params.input) as List<Sample>)
    RNASEQ(ch_samples)
}
```

Where `loadData()` is a generic data-loading function that supports multiple file formats (CSV, JSON, YAML, etc). The file contents must be compatible with the declared element type; an error is thrown if they are not. CSV files must include a header row and use a comma as the column separator. This function may be extended via plugin to support additional formats (e.g. Parquet).

The channel input can use a generic type such as `Map` or `Record`, or a custom record type to enable further validation. In the above example, using the `Sample` type ensures that each samplesheet row is validated against the record fields and the `fastq_1` and `fastq_2` columns are treated as file paths.

### Saving output channels to index files

The `emit:` section of a workflow can be treated like the `publish:` section of a entry workflow, defining which files are *terminal outputs* vs *intermediate outputs*. However, the output directory structure cannot be automatically inferred from the `emit:` section. It is normally specified by the output `path` directive, and does not necessarily correspond to the structure of the output channels.

When executing a named workflow directly, output files are not published to an output directory. Instead, the workflow output printed by Nextflow simply refers to output files by their work directory path.

This approach aligns with our goal to create a global content-addressable data store for files produced by Nextflow pipelines. A global data store with global search, caching, and automatic cleanup eliminates the need for per-run output directories. If needed, an output directory can be reconstructed from the workflow output printed by Nextflow.
