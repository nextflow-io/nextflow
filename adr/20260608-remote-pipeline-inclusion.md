# Remote pipeline inclusion

- Authors: Ben Sherman
- Status: draft
- Date: 2026-06-08
- Tags: pipelines, modules, dsl, registry

## Summary

Add the ability to include a remote pipeline into a *meta-pipeline*.

## Problem Statement

Nextflow supports reusing process definitions via remote *module* inclusion (e.g. `include { BWA_MEM } from 'nf-core/bwa/mem'`), but there is no standard mechanism to reuse an entire *pipeline* as a building block. Users must either fork and copy code, or compose/chain multiple `nextflow run` sessions which forfeits dataflow composition.

The natural unit of reuse for a pipeline is its *core workflow* -- the named workflow that takes and emits channels (e.g. `NFCORE_RNASEQ` in `nf-core/rnaseq`) -- as distinct from the deployment shell around it (the `params`, entry workflow, `output` block, and config). The nf-core community has already structured their pipelines around this split in anticipation of meta-pipelines.

The decision to make: how should a remote pipeline be distributed, resolved, stored, and included so that it can be composed into a meta-pipeline while preserving dataflow composition.

## Goals

- **Preserve dataflow composition**: the included pipeline participates in the meta-pipeline's dataflow graph (same session, same DAG, same work dir), enabling incremental reaction to emitted outputs.

- **Preserve reproducibility**: an included pipeline should produce the exact same results as it would when executed directly. Transitive dependencies should not be silently altered to reduce duplication.

- **Reuse existing conventions**: follow the conventions established by the module system (e.g. include syntax) as much as possible rather than introducing parallel conventions.

## Non-goals

- **Nested pipeline execution**: avoid Nextflow-in-Nextflow execution, which forfeits dataflow composition.

- **Pipeline execution via registry**: out of scope for first iteration. Registry-based execution (e.g. `nextflow pipeline run nf-core/rnaseq@3.0.0`) may be investigated in the future.

## Decision

Allow remote pipelines to be included into meta-pipelines, using the same namespacing conventions and include syntax as modules. Store the included pipeline in the meta-pipeline repository under `workflows/<scope>/<name>/` with its own subdirectories for modules and subworkflows. The meta-pipeline owns all top-level concerns (entry workflow, params, outputs, config). Pipelines should be written with a self-contained core workflow to make importing as easy as possible.

## Core Capabilities

### Composition over orchestration

The included pipeline must be incorporated into the meta-pipeline's dataflow graph in order to maximize dataflow concurrency. Existing approaches like pipeline chaining and Nextflow-in-Nextflow impose a synchronization barrier between each pipeline run.

To this end, the included pipeline should be written in a way that separates the *core workflow* from the rest of the pipeline (entry workflow, params, publishing, config). Only the core workflow is included into the meta-pipeline; the rest is discarded.

### Meta-pipeline owns all top-level concerns

The meta-pipeline owns the entry workflow, `params` block, `output` block, and config. An included pipeline contributes none of these. This approach aligns with existing include semantics, which only supports composition of processes and named workflows.

As a result, if a user wants to preserve any top-level concerns from the included pipeline, they must be explicitly replicated in the meta-pipeline. For example, params exposed by the included pipeline must be replicated in the meta-pipeline params and passed to the included pipeline's core workflow.

### Best practices for included pipeline

To be importable in practice, a pipeline's core workflow (and its dependent modules/workflows) should be free of external context:

1. Avoid `params` usage outside the entry workflow -- pass values as explicit process/workflow inputs.
2. Avoid `publishDir` -- use the `output` block.
3. Avoid use of project-level assets (`projectDir`, `bin`, `lib`) within the core workflow. Module-level assets can be safely used through the module `resources/` bundle and `moduleDir`.
4. Declare software dependencies (`container`, `conda`) in the process definition, not in config.
5. Avoid default `ext` settings in config -- specify these defaults in the process definition or use explicit process inputs. Otherwise, any default `ext` settings must be replicated manually in the meta-pipeline.
6. Avoid plugin functions within the core workflow.

For process directives, it is helpful to distinguish *what* is computed vs *how* it is computed. Directives that affect the *what* (`container`, `ext` settings) should be owned by the process definition. Directives that affect the *how* (`cpus`, `memory`, `executor`, `queue`, `errorStrategy`) should be owned by the meta-pipeline.

These constraints are not absolute -- it is possible to import a pipeline that does not adhere to any of these rules. Following these constraints simply makes it easier to import a pipeline with minimal extra work (manual replication, cross-cutting concerns).

### Pipeline inclusion and storage

Modules and pipelines share the same include syntax and naming conventions:

```groovy
// module
include { BWA_MEM } from 'nf-core/bwa/mem'

// pipeline
include { NFCORE_RNASEQ } from 'nf-core/rnaseq'
```

Including a remote pipeline is equivalent to including the top-level `main.nf` of that pipeline; any named workflow defined there can be included by name. By convention, the main script defines only the entry workflow and the core workflow (e.g. `NFCORE_RNASEQ` in nf-core/rnaseq). From this point, the included workflow can be called like any other workflow.

When a pipeline is included, it is vendored into the meta-pipeline project under `workflows/<scope>/<name>/`. Included pipelines are isolated -- each included pipeline has its own `modules/` and `workflows/` directories. This way, two pipelines can use different versions of the same module without compromising reproducibility.

Included pipelines should be committed to the meta-pipeline repository. The pipeline should have a *pipeline spec* (`nextflow_spec.json`) which specifies the pipeline version, so that Nextflow can track local changes.

## Open Questions

### Sourcing from Nextflow registry vs Git repositories

Pipelines could be stored in the Nextflow registry (as a new artifact type) or fetched directly from Git repositories. The pipeline registry is a potential long-term goal with other use cases, but it likely introduces additional scope that is not strictly related to meta-pipelines. Using existing Git repositories would be an expedient solution for the first iteration.

### Pipeline CLI

Sourcing remote pipelines from the Nextflow registry implies a `nextflow pipeline` command group for publishing and installing pipelines, similar to `nextflow module`. Even with a Git-based approach, a CLI is likely still needed to install and update remote pipelines.

### Using plugin functions in included pipeline

If an included pipeline uses plugin functions in the core workflow, these plugins must be explicitly declared in the meta-pipeline config, since the included pipeline config is not inherited.

Alternatively, these core plugin dependencies could be specified in the pipeline spec under `requires.plugins`. When installing a pipeline, Nextflow could copy these plugin declarations into the meta-pipeline config and/or spec.

Since this use case is rare -- plugin functions are typically used in the entry workflow outside the core workflow -- it can be deferred in the first iteration.

## Alternatives

### Pipeline chaining

An alternative to a meta-pipeline is a *pipeline chain*, in which multiple Nextflow pipelines are called in sequence via `nextflow run`.

For example, a fetchngs -> rnaseq pipeline chain can be implemented in a shell script:

```bash
# fetch FASTQ samples from NCBI SRA
nextflow -q run nf-core/fetchngs \
    --input samplesheet.csv \
    -output-format json \
    > results/output-fetchngs.json

# adapt fetchngs output to rnaseq input (add strandedness column)
nextflow -q run ./fetchngs-rnaseq.nf \
    -params-file results/output-fetchngs.json \
    --strandedness auto \
    -output-format json \
    > results/output-fetchngs-rnaseq.json

# perform RNAseq analysis
nextflow -q run nf-core/rnaseq \
    -params-file results/output-fetchngs-rnaseq.json \
    -output-format json \
    > results/output-rnaseq.json
```

Or a Nextflow pipeline:

```groovy
include { NEXTFLOW_RUN as NFCORE_FETCHNGS } from "./modules/local/nextflow/run"
include { NEXTFLOW_RUN as NFCORE_RNASEQ } from "./modules/local/nextflow/run"

workflow {
    NFCORE_FETCHNGS (
        'nf-core/fetchngs',
        // nextflow opts, pipeline inputs, etc ...
    )
    NFCORE_RNASEQ (
        'nf-core/rnaseq',
        // nextflow opts, pipeline inputs, etc ...
    )
}
```

The `NEXTFLOW_RUN` process simply calls `nextflow run` in a native process. See [nf-cascade](https://github.com/mahesh-panchal/nf-cascade) for more information about this approach.

Pipeline chains can also be implemented in Seqera Platform using actions (e.g. when a fetchngs run completes -> launch rnaseq on the fetchngs output).

Pipeline chaining works with any Nextflow pipeline out of the box, because it simply executes each pipeline directly. Chaining often requires glue logic to adapt upstream outputs to downstream outputs -- missing columns, different column names, etc -- but language features such as [workflow outputs](20251020-workflow-outputs.md) and [record types](20260306-record-types.md) make it easier by allowing pipelines to defined structured inputs and outputs.

The downside of pipeline chaining is that it sacrifices dataflow concurrency -- each pipeline must complete before the next pipeline can start. As a result, pipeline chaining is a categorically different solution from meta-pipelines. It remains a valid option for certain use cases, such as simple chains (A -> B -> C) of off-the-shelf pipelines. For compositions that are more complex and/or require maximum dataflow concurrency, meta-pipelines are the general solution.

### Best of both: runtime inheritance

The meta-pipeline approach treats the included pipeline as a *white box* -- it composes the pipeline like any included workflow, producing a single dataflow graph. It also imposes several constraints on how the included pipeline is written, and it imposes development overhead. For example, any params / outputs that need to be exposed from the included pipeline must be replicated in the meta-pipeline.

Pipeline chaining treats the included pipeline as a *black box* -- it preserves the exact pipeline behavior (core workflow + entry workflow + config) while forfeiting dataflow composition (separate dataflow graphs).

An ideal solution might combine the best of both: compose pipelines into a single dataflow graph (white box) while inheriting each pipeline's params, outputs, and config so they need not be replicated (black box). We considered such a model, where an included pipeline contributes its shell as namespaced, overridable defaults, but rejected it. Dataflow composition fundamentally requires exposing the core workflow as a set of channel ports, so the white-box mechanism is unavoidable; inheritance would only layer implicit behavior on top of it. That behavior comes at a steep cost: it relocates a one-time *write* cost (boilerplate) into a recurring *read* cost (hidden defaults, auto-bound arguments, auto-published outputs), burdens every tool that must now understand it (linter, type checker, config resolution, resume), and conflicts with the frozen-island philosophy that otherwise governs vendored code.

Instead, we keep the meta-pipeline fully explicit and address the boilerplate at write time. The replicated params, outputs, and workflow call are mechanical transcriptions of the included pipeline's shell -- precisely the kind of task a coding agent can generate from the pipeline definition and a description of which params and outputs to expose, leaving the developer to write only the composition logic that carries novel intent.

## Links

- Related: [Module system](20251114-module-system.md)
- Related: [Workflow params](20250825-workflow-params.md)
- Related: [Workflow outputs](20251020-workflow-outputs.md)

## Appendix

### Example: fetchngs -> rnaseq

This section walks through the aforementioned `fetchngs -> rnaseq` example as a meta-pipeline.

> NOTE: This example uses simplified and idealized versions of `nf-core/fetchngs` and `nf-core/rnaseq` and may not match the actual implementations.

**Project layout**

The meta-pipeline is an ordinary Nextflow project with `nf-core/fetchngs` and `nf-core/rnaseq` vendored under `workflows/`:

```
fetchngs-rnaseq/
├── main.nf
├── nextflow.config
└── workflows/
    └── nf-core/
        ├── fetchngs/
        │   ├── main.nf
        │   ├── nextflow_spec.json
        │   ├── modules/
        │   └── workflows/
        └── rnaseq/
            ├── main.nf
            ├── nextflow_spec.json
            ├── modules/
            └── workflows/
```

Each pipeline has its own `modules/` and `workflows/`, so the two pipelines can depend on different versions of the same module without conflict. Both pipelines are committed to the meta-pipeline repository.

**Pipeline code**

The included pipelines are defined as follows, with a clear separation of *core workflow* from *entry workflow*:

```groovy
// nf-core/fetchngs — main.nf
params {
    input: Path // file of SRA/ENA accessions
}
workflow {
    main:
    ch_ids = channel.fromPath(params.input).splitCsv()
    ch_samples = NFCORE_FETCHNGS( ch_ids )
    publish:
    samples = ch_samples
}
output {
    samples: Channel<Sample> { path 'fastq' }
}

workflow NFCORE_FETCHNGS {
    take:
    ids: Channel<String>

    main:
    // ...

    emit:
    samples: Channel<Sample>
}
```

```groovy
// nf-core/rnaseq — main.nf
params {
    input: Path // samplesheet
    aligner: String = 'star_salmon'
    fasta: Path
}
workflow {
    main:
    ch_samples = channel.fromPath(params.input).splitCsv()
    rnaseq = NFCORE_RNASEQ( ch_samples, params.aligner, params.fasta )
    publish:
    multiqc = rnaseq.multiqc
    bams    = rnaseq.bams
    counts  = rnaseq.counts
}
output {
    multiqc: Path { path 'multiqc' }
    bams: Channel<Path> { path 'bams' }
    counts: Channel<Path> { path 'counts' }
}

workflow NFCORE_RNASEQ {
    take:
    samples: Channel<Sample>
    aligner: String
    fasta: Path

    main:
    // ...

    emit:
    multiqc: Value<Path>
    bams: Channel<Path>
    counts: Channel<Path>
}
```

The meta-pipeline includes the core workflow from each pipeline and composes them into an entry workflow with params and outputs:

```groovy
include { NFCORE_FETCHNGS } from 'nf-core/fetchngs'
include { NFCORE_RNASEQ } from 'nf-core/rnaseq'

params {
    input: Path
    strandedness: String = 'auto'
    aligner: String = 'star_salmon'
    fasta: Path
}

workflow {
    main:
    // fetch FASTQ samples from NCBI SRA
    ch_ids = channel.fromPath(params.input).splitCsv()
    ch_samples = NFCORE_FETCHNGS( ch_ids )

    // adapt fetchngs output to rnaseq input (add strandedness)
    ch_samples = ch_samples.map { r ->
        r + record(strandedness: params.strandedness)
    }

    // perform RNAseq analysis
    rnaseq = NFCORE_RNASEQ( ch_samples, params.aligner, params.fasta )

    publish:
    multiqc = rnaseq.multiqc
    bams    = rnaseq.bams
    counts  = rnaseq.counts
}

output {
    multiqc: Path { path 'multiqc' }
    bams: Channel<Path> { path 'bams' }
    counts: Channel<Path> { path 'counts' }
}
```

Notes about the white-box approach:

- **The handoff is a channel, not a file.** A pipeline chain blocks until fetchngs finishes before rnaseq starts. Here, `ch_samples` is a live channel: rnaseq begins aligning each sample the moment fetchngs emits it. This is the dataflow composition that motivates the meta-pipeline.

- **The adapter is an operator, not a pipeline.** The strandedness gap that required a separate `fetchngs-rnaseq.nf` adapter in the chaining example collapses to a single `map` operator in the meta-pipeline.

- **Params and outputs are replicated, not inherited.** `--input` and `--strandedness` are declared in the meta-pipeline's own `params` block and passed explicitly into the core workflows. Similarly, any outputs must be declared as such in the meta-pipeline's `output` block. The included pipelines do not contribute any of their own params, entry workflows, or output blocks.

**Configuration**

Since each included pipeline is just part of the dataflow graph, configuration works like normal. Processes in an included pipeline can be targeted via config selector:

```groovy
process {
    withName: 'NFCORE_FETCHNGS:.*:SRATOOLS_FASTERQDUMP' {
        cpus   = 6
        memory = 24.GB
    }
    withName: 'NFCORE_RNASEQ:.*:STAR_ALIGN' {
        cpus   = 12
        memory = 72.GB
    }
}
```

Both the meta-pipeline developer and users can override whatever they want from config. In practice, the process definitions should own the *what* (`container`, `conda`) while the meta-pipeline config should own the *how* (`cpus`, `memory`).

**Trade-offs**

| Concern | Chaining (black-box) | Meta-pipeline (white-box) |
| --- | --- | --- |
| Concurrency | Synchronous (each `run` completes first) | Asynchronous (rnaseq reacts to each fetchngs sample) |
| Params and outputs | Owned by each pipeline | Replicated in the meta-pipeline |
| Resource config | Per-pipeline config files | Unified meta-pipeline config |

The most notable trade-off is the replication of params and outputs: anything the included pipelines exposed at the top level (params, published outputs) must be re-declared in the meta-pipeline.
