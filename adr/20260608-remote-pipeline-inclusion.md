# Remote pipeline inclusion

- Authors: Ben Sherman
- Status: draft
- Date: 2026-06-08
- Tags: pipelines, modules, dsl, registry
- Version: 1.1

## Updates

### Version 1.1 (2026-06-22)
- **Separate remote pipelines from remote workflows**: Workflows are treated separately by the [Workflow modules ADR](./20260608-workflow-modules.md).
- **Replace core workflow distinction with pipeline inclusion**: Instead of isolating the *core workflow* of a pipeline, the include syntax is extended to support *pipeline inclusion*, in which the `params` / `workflow` / `output` trio is imported and used like a named workflow.

## Summary

Add the ability to include a remote pipeline into a *meta-pipeline*.

## Problem Statement

Nextflow supports reusing process definitions via remote *module* inclusion (e.g. `include { BWA_MEM } from 'nf-core/bwa/mem'`), but there is no standard mechanism to reuse an entire *pipeline* as a building block. Users must either fork and copy code, or compose/chain multiple `nextflow run` sessions which forfeits dataflow composition.

The [module system](20251114-module-system.md) and [workflow module](20260608-workflow-modules.md) ADRs define how standalone *processes* and *workflows* should be distributed as modules through the Nextflow registry. This ADR defines how *pipelines* -- workflows with a deployment shell -- should be composed into larger *meta-pipelines*.

## Goals

- **Preserve dataflow composition**: the included pipeline participates in the meta-pipeline's dataflow graph (same session, same DAG, same work dir), enabling incremental reaction to emitted outputs.

- **Preserve reproducibility**: an included pipeline should produce the exact same results as it would when executed directly. Transitive dependencies should not be silently altered to reduce duplication.

- **Reuse existing conventions**: follow the conventions established by the module system (e.g. include syntax) as much as possible rather than introducing parallel conventions.

## Non-goals

- **Nested pipeline execution**: avoid Nextflow-in-Nextflow execution, which forfeits dataflow composition.

- **Pipeline execution via registry**: out of scope for first iteration. Registry-based execution (e.g. `nextflow pipeline run nf-core/rnaseq@3.0.0`) may be investigated in the future.

## Decision

Allow pipelines to be published and installed through the Nextflow registry, using the same namespacing conventions as modules. Store the included pipeline in the meta-pipeline repository under `pipelines/<scope>/<name>/` with its own subdirectories for modules and workflows. Provide a way to include an entire pipeline (`params` block, entry workflow, `output` block) as a named workflow to facilitate workflow composition.

## Core Capabilities

### Pipeline composition

A pipeline -- that is, a `params` / `workflow` / `output` trio -- can be included and called like a named workflow. This way, pipelines can be composed using regular dataflow logic.

For example, given the following pipeline:

```groovy
// pipelines/rnaseq.nf
params {
    input: Path
    aligner: String = 'star_salmon'
    fasta: Path
}
workflow {
    // ...
}
output {
    bams: Channel<Path> { path 'bams' }
    multiqc: Path { path 'multiqc' }
}
```

It can be included and called as follows:

```groovy
// main.nf
include { workflow as RNASEQ } from './pipelines/rnaseq.nf'

workflow {
    rnaseq = RNASEQ(
        input: file('input.csv'),
        fasta: file('index.fasta')
    )
    rnaseq.bams.view()      // Channel<Path>
    rnaseq.multiqc.view()   // Value<Path>
}
```

Notes:

- The pipeline must included using the `workflow` keyword and aliased to a specific name (`RNASEQ`).
- The `params` block becomes the `take:` section and the `output` block becomes the `emit:` section.
- The workflow is called using named arguments so that defaults can be omitted.
- All outputs are either a `Channel` or wrapped as `Value<T>`, allowing them to be used in regular dataflow logic.

### Remote pipeline inclusion and storage

Pipelines can be published, installed, and included through the Nextflow registry:

```groovy
// module
include { BWA_MEM } from 'nf-core/bwa/mem'

// pipeline
include { workflow as NFCORE_RNASEQ } from 'nf-core/rnaseq'
```

When a pipeline is included from the registry, it is vendored into the including project under `pipelines/<scope>/<name>/`. Included pipelines are isolated -- each included pipeline has its own `modules/` directory. This way, two pipelines can use different versions of the same module without compromising reproducibility.

Included pipelines should be committed to the meta-pipeline repository. The pipeline should have a *pipeline spec* (`nextflow_spec.json`) which specifies the pipeline version, so that Nextflow can track local changes.

### Best practices for including pipelines

Pipeline inclusion only captures the pipeline's main script and included modules -- it does not capture external context such as config or the `lib` directory. As a result, the pipeline should be written in a way that works when included by a meta-pipeline:

1. Pipeline parameters should be defined in the script `params` block. The config should only declare *config params* (params that only affect config settings).

2. Project-level assets (`projectDir`, `bin`, `lib`) should not be used since the meta-pipeline will have a different project root. Module-level assets can be safely used through the module `resources/` bundle and `moduleDir`.

3. Default `ext` settings should be specified in the process definition or avoided in favor of process inputs.

4. Software dependencies (`container`, `conda`) should be declared in the process definition, not in config.

5. Workflow outputs should be published using the `output` block, not `publishDir`.

None of these constraints are absolute. All of them can be circumvented by manually replicating the external context in the meta-pipeline. Following these constraints simply makes it easier to import a pipeline with minimal extra work.

## Open Questions

### Pipeline registry and CLI

Sourcing remote pipelines from the Nextflow registry implies a pipeline registry API and a `nextflow pipeline` command group for publishing and installing pipelines. This infrastructure can be largely inferred from existing patterns established for modules.

One aspect that remains open is the pipeline spec (`nextflow_spec.json` or `nextflow_schema.json`) which may have a different shape from the module spec (`meta.yml`). A minimal pipeline spec could be introduced to enable remote pipelines without bloating scope. Alternatively, the pipeline version could be managed by a helper file (e.g. `.pipeline-info`) until the pipeline spec is finalized.

### Using plugin functions in included pipeline

If an included pipeline uses plugins, these plugins must be explicitly declared in the meta-pipeline config since they cannot be inferred from the pipeline inclusion.

Alternatively, these core plugin dependencies could be specified in the pipeline spec under `requires.plugins`. When installing a pipeline, Nextflow could copy these plugin declarations into the meta-pipeline config and/or spec.

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

params {
    // ...
}

workflow {
    // fetch FASTQ samples from NCBI SRA
    fetchngs = NFCORE_FETCHNGS (
        'nf-core/fetchngs',
        // nextflow opts, pipeline inputs, etc ...
    )
    // adapt fetchngs output to rnaseq input (add strandedness column)
    ch_samples = fetchngs2rnaseq(fetchngs)
    // perform RNAseq analysis
    rnaseq = NFCORE_RNASEQ (
        'nf-core/rnaseq',
        // nextflow opts, pipeline inputs, etc ...
    )
}

output {
    // ...
}
```

The `NEXTFLOW_RUN` process simply calls `nextflow run` in a native process. See [nf-cascade](https://github.com/mahesh-panchal/nf-cascade) for more information about this approach.

Pipeline chains can also be implemented in Seqera Platform using actions (e.g. when a fetchngs run completes -> launch rnaseq on the fetchngs output).

Pipeline chaining works with any Nextflow pipeline out of the box, because it simply executes each pipeline directly. Language features such as [workflow outputs](20251020-workflow-outputs.md) and [record types](20260306-record-types.md) make pipeline chaining easier by allowing each pipeline to define structured inputs and outputs.

However, there are a number of downsides:

- It forfeits native dataflow composition. The developer must serialize/deserialize samplesheet files instead of passing channels directly between pipelines. Each pipeline must complete before the next pipeline can start.

- It requires an external workflow system instead of reusing the language that pipeline developers already know. Even the Nextflow-in-Nextflow approach shown above requires many tricks to orchestrate nested pipeline runs via the `NEXTFLOW_RUN` process.

Pipeline chaining can be practical for certain use cases, such as simple chains (A -> B -> C) of off-the-shelf pipelines. But the general solution is to compose pipelines using dataflow logic, just like any other Nextflow pipeline.

## Links

- Related: [Module system](20251114-module-system.md)
- Related: [Workflow params](20250825-workflow-params.md)
- Related: [Workflow outputs](20251020-workflow-outputs.md)

## Appendix

### Example: fetchngs -> rnaseq

This section walks through the aforementioned `fetchngs -> rnaseq` example as a meta-pipeline.

> NOTE: This example uses simplified and idealized versions of `nf-core/fetchngs` and `nf-core/rnaseq` and may not match the actual implementations.

**Project layout**

The meta-pipeline is an ordinary Nextflow project with `nf-core/fetchngs` and `nf-core/rnaseq` vendored under `pipelines/`:

```
fetchngs-rnaseq/
├── main.nf
├── nextflow.config
└── pipelines/
    └── nf-core/
        ├── fetchngs/
        │   ├── main.nf
        │   ├── nextflow_spec.json
        │   └── modules/
        └── rnaseq/
            ├── main.nf
            ├── nextflow_spec.json
            └── modules/
```

Each pipeline has its own `modules/` directory, so the two pipelines can depend on different versions of the same module without conflict. Both pipelines are committed to the meta-pipeline repository.

**Pipeline code**

The included pipelines are defined as follows:

```groovy
// nf-core/fetchngs — main.nf
params {
    input: Path // file of SRA/ENA accessions
}
workflow {
    main:
    ch_ids = channel.fromPath(params.input).splitCsv()
    ch_samples = // ...
    publish:
    samples = ch_samples
}
output {
    samples: Channel<Sample> { path 'fastq' }
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
    rnaseq = // ...
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

The meta-pipeline includes each pipeline and composes them into a new entry workflow with params and outputs:

```groovy
include { workflow as NFCORE_FETCHNGS } from 'nf-core/fetchngs'
include { workflow as NFCORE_RNASEQ } from 'nf-core/rnaseq'

params {
    input: Path
    strandedness: String = 'auto'
    aligner: String = 'star_salmon'
    fasta: Path
}

workflow {
    main:
    // fetch FASTQ samples from NCBI SRA
    fetchngs = NFCORE_FETCHNGS( input: params.input )

    // adapt fetchngs output to rnaseq input (add strandedness)
    ch_samples = fetchngs.samples.map { r ->
        r + record(strandedness: params.strandedness)
    }

    // perform RNAseq analysis
    rnaseq = NFCORE_RNASEQ( input: ch_samples, aligner: params.aligner, fasta: params.fasta )

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

Notes:

- **The handoff is a channel, not a file.** A pipeline chain blocks until fetchngs finishes before rnaseq starts. Here, `ch_samples` is a live channel: rnaseq begins aligning each sample the moment fetchngs emits it. This is the dataflow composition that motivates the meta-pipeline.

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

Both the meta-pipeline developer and users can override whatever they want from config.
