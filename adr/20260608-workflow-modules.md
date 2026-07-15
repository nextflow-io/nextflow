# Workflow modules

- Authors: Ben Sherman
- Status: draft
- Date: 2026-06-08
- Tags: workflows, modules, registry

## Summary

Add the ability to include a workflow as a module from the Nextflow registry.

## Problem Statement

A Nextflow *module* is currently defined as a standalone process definition (with corresponding spec file). The [module system ADR](20251114-module-system.md) defines how these modules are published, distributed, and executed through the Nextflow registry.

There is a similar need to share and re-use workflows / subworkflows. The nf-core community has curated a collection of re-usable subworkflows in the [nf-core/modules](https://github.com/nf-core/modules) repository.

The module system should be extended to include both standalone *processes* and *workflows*.

## Goals

- **Workflow inclusion**: allow workflows to be included from the registry using the same include syntax as for modules.

- **Workflow execution**: allow workflows to be executed directly via `nextflow module run`.

- **Reuse existing infrastructure**: the module system already defines a registry API and CLI for publishing and installing modules. Extend this system to include workflows instead of treating workflows as a separate concept.

## Non-goals

- **Pipeline inclusion**: *pipelines* are distinct from *workflows* -- they specify params, publishing, and config, not just workflow logic.

## Decision

Extend the definition of *module* to include both standalone *processes* and standalone *workflows*. Allow workflow modules to be published, installed, included, and executed via the Nextflow registry. Store installed workflows under `modules/<scope>/<name>/` alongside process modules. Allow each workflow module to have its own `modules` directory for transitive dependencies.

## Core Capabilities

### Workflow modules vs process modules

A module can refer to a *workflow module* or *process module*, depending on whether it defines a workflow or process.

Workflow modules can be published to and queried from the Nextflow registry using the same modules API and `nextflow module` CLI. Workflow modules also have the same directory structure as process modules.

Since workflow modules are just modules, a process module and workflow module cannot have the same name in a module namespace.

### Module spec

The module spec is extended as follows in order to support standalone workflows:

- A `kind` field to distinguish between workflow modules (`kind: Workflow`) and process modules (`kind: Process`).

- A `requires.modules` field to specify *transitive dependencies*, since a workflow can depend on other processes and workflows.

- Fields that are relevant only to process modules (i.e. `topics`, `tools`) cannot be specified for workflow modules.

For example, [nf-core/fastq_align_star](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/fastq_align_star):

```yaml
name: nf-core/fastq_align_star
kind: Workflow
version: 0.0.0-4e3e10e
description: Align reads to a reference genome using STAR then sort, index, and compute statistics with samtools
keywords:
  - align
  - fasta
  - genome
  - reference
authors:
  - "@JoseEspinosa"
license: MIT
requires:
  nextflow: ">=24.04.0"
  modules:
    - nf-core/star/align@0.0.0-4e3e10e
    - nf-core/bam_sort_stats_samtools@0.0.0-4e3e10e
```

Inputs and outputs are derived from the workflow's `take:` and `emit:` sections. For example (outputs truncated for brevity):

```yaml
input:
  - name: ch_reads
    type: channel
    description: List of input FastQ files of size 1 and 2 for single-end and paired-end data, respectively.
  - name: index
    type: directory
    description: STAR genome index
  - name: gtf
    type: file
    description: GTF file used to set the splice junctions with the --sjdbGTFfile flag
  - name: star_ignore_sjdbgtf
    type: boolean
    description: If true the --sjdbGTFfile flag is set
  - name: fasta_fai
    type: file
    description: Reference genome fasta file and index
  - name: transcripts_fasta_fai
    type: file
    description: Optional reference genome fasta file and index
output:
  - name: bam
    type: channel
    description: BAM file ordered by samtools
  - name: bai
    type: channel
    description: BAI index of the ordered BAM file
  - name: stats
    type: channel
    description: File containing samtools stats output
  - name: flagstat
    type: channel
    description: File containing samtools flagstat output
  - name: idxstats
    type: channel
    description: File containing samtools idxstats output
```

### Workflow inclusion and storage

Workflow modules use the same include syntax and naming conventions as process modules:

```groovy
// process
include { BWA_MEM } from 'nf-core/bwa/mem'

// workflow
include { FASTQ_ALIGN_STAR } from 'nf-core/fastq_align_star'
```

When a workflow is included, it is vendored into the including project under `modules/<scope>/<name>/`. Included workflows should be committed to the including repository.

Each workflow module should store its dependencies within its own `modules` directory. This way, two workflows can use different versions of the same module without introducing a version conflict. Dependencies are not bundled with the workflow module -- they are installed when the workflow module is installed.

For example, installing `nf-core/fastq_align_star` would produce the following directory tree:

```
my-pipeline/
├── main.nf                                 # include { FASTQ_ALIGN_STAR } from 'nf-core/fastq_align_star'
├── nextflow.config
└── modules/
    └── nf-core/
        └── fastq_align_star/
            ├── .module-info
            ├── main.nf                      # includes STAR_ALIGN, BAM_SORT_STATS_SAMTOOLS
            ├── meta.yml
            └── modules/
                └── nf-core/
                    ├── star/align/
                    │   ├── .module-info
                    │   ├── main.nf
                    │   └── meta.yml
                    └── bam_sort_stats_samtools/
                        ├── .module-info
                        ├── main.nf          # includes SAMTOOLS_SORT, SAMTOOLS_INDEX, SAMTOOLS_STATS, ...
                        ├── meta.yml
                        └── modules/
                            └── nf-core/
                                ├── samtools/sort/
                                │   ├── .module-info
                                │   ├── main.nf
                                │   └── meta.yml
                                ├── samtools/index/
                                │   ├── .module-info
                                │   ├── main.nf
                                │   └── meta.yml
                                ├── samtools/stats/
                                │   ├── .module-info
                                │   ├── main.nf
                                │   └── meta.yml
                                ├── samtools/idxstats/
                                │   ├── .module-info
                                │   ├── main.nf
                                │   └── meta.yml
                                └── samtools/flagstat/
                                    ├── .module-info
                                    ├── main.nf
                                    └── meta.yml
```

Note that `fastq_align_star` declares only its *direct* dependencies (`star/align` and `bam_sort_stats_samtools`) in `requires.modules`. The `samtools` processes are transitive dependencies of `bam_sort_stats_samtools` and are vendored within *its* own `modules` directory, not directly under `fastq_align_star`. Because each workflow module vendors its own dependencies, the same module may appear more than once in the tree when it is used by multiple workflows.

### Workflow execution

Typed workflows can be executed directly by inferring the `params` and `output` blocks from the `take:` and `emit:` sections.

For example, given the following workflow:

```groovy
workflow RNASEQ {
    take:
    samples: Channel<Sample>
    index: Path

    main:
    ch_aligned = ALIGN(samples, index)
    multiqc_report = MULTIQC(ch_aligned.collect())

    emit:
    aligned: Channel<AlignedSample> = ch_aligned
    multiqc_report: Path = multiqc_report
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}

record AlignedSample {
    id: String
    bam: Path
    bai: Path
}
```

The user can run the workflow directly as follows:

```bash
nextflow module run rnaseq.nf \
    --samples input.csv \
    --index index.fasta
```

Nextflow executes the `RNASEQ` workflow as if it were wrapped in the following entry workflow:

```groovy
params {
    samples: Channel<Sample>
    index: Path
}

workflow {
    main:
    rnaseq = RNASEQ(params.samples, params.index)

    publish:
    aligned = rnaseq.aligned
    multiqc_report = rnaseq.multiqc_report
}

output {
    aligned: Channel<AlignedSample> {}
    multiqc_report: Path {}
}
```

This way, the user can run a workflow directly without having to write an entry workflow for it.

Direct execution requires the ability to load an input channel directly from an index file (samplesheet). This can be done by loading the samplesheet data based on the file extension (CSV, JSON, YAML), casting each record to the given record type, and loading the collection as a channel. The data-loading function may be extended via plugin to support additional formats (e.g. Parquet).

The channel input can use a generic type such as `Map` or `Record`, or a custom record type to enable further validation. In the above example, using the `Sample` type ensures that each samplesheet row is validated against the record fields and the `fastq_1` and `fastq_2` columns are treated as file paths.

When executing a named workflow directly, output files are not published to an output directory. Instead, the workflow output printed by Nextflow simply refers to output files by their work directory path.

## Alternatives

### Workflows vs processes

One alternative is to treat workflows as a separate concept from modules, restricting the definition of *module* to only include standalone processes.

However, the broader meaning of *module* is a re-usable component, and both processes and workflows are re-usable components, so it makes more sense to extend the module system rather than introduce a parallel system for workflows. A parallel system would also require duplicating a lot of existing code (registry API, CLI, etc).

Process modules and workflow modules can be distinguished by different specializations of the module spec (`kind: Process` vs `kind: Workflow`).

### Workflows vs subworkflows

The nf-core community makes a distinction between *workflows* and *subworkflows*:

- Subworkflows are stored in the `subworkflows/` directory and can be shared across pipelines.

- Workflows are stored in the `workflows/` directory and are owned by a specific pipeline.

Nextflow makes no such distinction, as there is no functional difference between a workflow and a subworkflow. Both are stored in the `modules/` directory and both can be published and installed through the Nextflow registry.

The Nextflow registry makes it easier for pipelines to share and reuse workflows, since the pipeline can be the source-of-truth rather than having to surrender ownership to the `nf-core/modules` repository.

### Workflows vs pipelines

One alternative is to treat all workflows as pipelines, restricting the definition of *module* to only include standalone processes.

Aside from the problems mentioned under [Workflows vs processes](#workflows-vs-processes), this approach conflates two commonly-understood concepts: a *pipeline* which is an end-to-end analysis, and a *workflow* or *subworkflow* which is a component in a larger analysis. It is more congruent with conventional understanding to treat workflows (subworkflows) as a separate concept from pipelines.

### Transitive dependencies

Workflow modules, unlike process modules, can include other modules (both processes and workflows). As a result, we must address the *diamond dependency problem* -- what happens when two workflows depend on different versions of the same module?

Currently, nf-core enforces a flat module structure. nf-core modules are currently versioned by commit hash and workflows pin exact versions. There is no way to resolve conflicts between workflows. Instead, all nf-core workflows are kept in sync at publish time. This approach places the maintenance burden on *module developers* -- when updating a module, you must also update all consuming workflows.

Instead, the ADR avoids this problem entirely by allowing each workflow to have its own `modules` directory. The downside is that modules will be duplicated in the project, making the project larger and increasing the number of scripts to be parsed.

## Links

- [nf-core terminology](https://nf-co.re/docs/usage/getting_started/terminology)
- Community discussions: [#4112](https://github.com/nextflow-io/nextflow/issues/4112), [nf-core/modules#8](https://github.com/nf-core/modules/issues/8)
- Refines: [Module system](20251114-module-system.md)
