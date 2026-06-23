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

- **Reuse existing infrastructure**: the module system already defines a registry API and CLI for publishing and installing modules. Treating workflows as a separate concept would require duplicating much of this infrastructure.

## Non-goals

- **Pipeline inclusion**: *pipelines* are distinct from *workflows* -- they specify params, publishing, and config, not just workflow logic.

## Decision

Extend the definition of *module* to include both standalone *processes* and standalone *workflows*. Allow workflow modules to be published, installed, included, and executed via the Nextflow registry. Store the included workflow in the including repository under `workflows/<scope>/<name>/`.

## Core Capabilities

### Workflow modules vs process modules

A module can refer to a *workflow module* or *process module*, depending on whether it defines a workflow or process.

Workflow modules can be published to and queried from the Nextflow registry using the same modules API and `nextflow module` CLI. Workflow modules also have the same directory structure as process modules.

Since workflow modules are just modules, a process module and workflow module cannot have the same name in a module namespace.

### Module spec

The module spec is extended as follows in order to support standalone workflows:

- A `kind` field to distinguish between workflow modules (`kind: Workflow`) and process modules (`kind: Process`).

- A `requires.modules` field to specify *transitive dependencies*, since a workflow can depend on other processes and workflows.

- Certain fields (`topics`, `tools`) cannot be specified for workflow modules because they are not applicable.

For example:

```yaml
name: nf-core/fastq_align_star
kind: Workflow
version: 0.0.0-4e3e10e
description: Align reads to a reference genome using bowtie2 then sort with samtools
authors:
  - "@JoseEspinosa"
license: MIT
requires:
  nextflow: ">=24.04.0"
  modules:
    - nf-core/star/align@0.0.0-4e3e10e
    - nf-core/samtools/sort@0.0.0-4e3e10e
    - nf-core/samtools/index@0.0.0-4e3e10e
    - nf-core/samtools/stats@0.0.0-4e3e10e
    - nf-core/samtools/idxstats@0.0.0-4e3e10e
    - nf-core/samtools/flagstat@0.0.0-4e3e10e
    - nf-core/bam_sort_stats_samtools@0.0.0-4e3e10e
```

### Workflow inclusion and storage

Workflow modules use the same include syntax and naming conventions as process modules:

```groovy
// process
include { BWA_MEM } from 'nf-core/bwa/mem'

// workflow
include { FASTQ_ALIGN_STAR } from 'nf-core/fastq_align_star'
```

When a workflow is included, it is vendored into the including project under `workflows/<scope>/<name>/`. Included workflows should be committed to the including repository.

Transitive dependencies (specified by `requires.modules`) should also be installed in the `modules/` and `workflows/` directories alongside the included workflow. Since modules are flattened, it is not possible for a pipeline to use two different versions of the same process or workflow.

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

Process modules and workflow modules can be distinguished by different specializations of the module spec (`kind: Process` vs `kind: Workflow`) and different storage locations (`modules/` vs `workflows/`).

### Workflows vs subworkflows

The nf-core community makes a distinction between *workflows* and *subworkflows*:

- Subworkflows are stored in the `subworkflows/` directory and can be shared across pipelines.

- Workflows are stored in the `workflows/` directory and are owned by a specific pipeline.

Nextflow makes no such distinction, as there is no functional difference between a workflow and a subworkflow. Both are stored in the `workflows/` directory and both can be published and installed through the Nextflow registry.

The Nextflow registry makes it easier for pipelines to share and reuse workflows, since the pipeline can be the source-of-truth rather than having to surrender ownership to the `nf-core/modules` repository.

### Workflows vs pipelines

One alternative is to treat all workflows as pipelines, restricting the definition of *module* to only include standalone processes.

Aside from the problems mentioned under [Workflows vs processes](#workflows-vs-processes), this approach conflates two commonly-understood concepts: a *pipeline* which is an end-to-end analysis, and a *workflow* or *subworkflow* which is a component in a larger analysis. It is more congruent with conventional understanding to treat workflows (subworkflows) as a separate concept from pipelines.

## Links

- [nf-core terminology](https://nf-co.re/docs/usage/getting_started/terminology)
- Refines: [Module system](20251114-module-system.md)
- Related: [Remote pipeline inclusion](20260608-remote-pipeline-inclusion.md)
