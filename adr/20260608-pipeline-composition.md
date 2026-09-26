# Pipeline composition

- Authors: Ben Sherman
- Status: proposed
- Date: 2026-06-08
- Tags: pipelines, modules, dsl
- Version: 1.4

## Updates

### Version 1.4 (2026-09-25)

- **Remove output block inclusion**: the `output` block of a pipeline can no longer be included. A pipeline call returns its outputs like a workflow call, and the meta-pipeline declares its outputs like any other pipeline.

### Version 1.3 (2026-09-17)

- **Included output block is only a record type**: including the `output` block of a pipeline provides a record type of its outputs and nothing else, mirroring the `params` block. Declaring an output with this type no longer redeclares each output of the included pipeline with its output directives. The meta-pipeline declares its outputs like any other pipeline.

### Version 1.2 (2026-07-13)

- **Reframe as pipeline composition**: the core feature is the ability to compose pipelines in a Nextflow-native manner. Meta-pipelines are the artifact. Remote pipeline inclusion is deferred to future work.

### Version 1.1 (2026-06-22)

- **Separate remote pipelines from remote workflows**: Workflows are treated separately by the [Workflow modules ADR](20260608-workflow-modules.md).
- **Replace core workflow distinction with pipeline inclusion**: Instead of isolating the *core workflow* of a pipeline, the include syntax is extended to support *pipeline inclusion*, in which the `params` / `workflow` / `output` trio is imported and used like a named workflow.

## Summary

Provide a way to compose pipelines using regular dataflow logic.

## Problem Statement

Processes and workflows can be composed into a larger workflow using dataflow logic. However, pipelines cannot be composed in the same way. The only way to call a pipeline is via `nextflow run`, which does not allow for dataflow composition.

This ADR defines how a pipeline can be included like a named workflow and composed with other pipelines with dataflow logic.

## Goals

- **Preserve dataflow composition**: the included pipeline participates in the including pipeline's dataflow graph (same session, same DAG, same work dir), enabling incremental reaction to emitted outputs.

- **Preserve reproducibility**: an included pipeline should produce the exact same results as it would when executed directly. Transitive dependencies should not be silently altered to reduce duplication.

## Non-goals

- **Nested pipeline execution**: avoid Nextflow-in-Nextflow execution, which forfeits dataflow composition.

- **Remote pipeline inclusion**: deferred to future work.

- **Remote pipeline execution**: out of scope. Registry-based execution (e.g. `nextflow pipeline run nf-core/rnaseq@3.0.0`) may be investigated in the future.

## Decision

Provide a way to include an entire pipeline (`params` block, entry workflow, `output` block) as a named workflow to facilitate workflow composition.

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
    rnaseq = RNASEQ(record(
        input: file('input.csv'),
        fasta: file('index.fasta')
    ))
    rnaseq.bams.view()      // Channel<Path>
    rnaseq.multiqc.view()   // Value<Path>
}
```

Notes:

- The pipeline must be included using the `workflow` keyword and aliased to a specific name (`RNASEQ`).
- A pipeline with a `params` block is called with a single record of params, so that params with a default can be omitted. A pipeline without a `params` block is called with no arguments.
- The `output` block becomes the `emit:` section. As with a workflow call, a pipeline with a single output returns it directly, and a pipeline with multiple outputs returns a record.
- All outputs are either a `Channel` or wrapped as `Value<T>`, allowing them to be used in regular dataflow logic.

### Including the params block

The `params` block of a pipeline can be included as a *record type*, so that a calling pipeline can declare the params of an included pipeline as a single param instead of redeclaring each one:

```groovy
// main.nf
include {
    params as RnaseqParams ;
    workflow as RNASEQ
} from './pipelines/rnaseq.nf'

params {
    rnaseq: RnaseqParams
}

workflow {
    rnaseq = RNASEQ( params.rnaseq )
}
```

Notes:

- The included params record type is *partial*: every field is nullable, and defaults are not pre-filled. The included pipeline applies its own defaults and validates its required params when it is called.
- The user can provide each param of the included pipeline as `--rnaseq.<name>`, and the calling pipeline can override specific params, e.g. `params.rnaseq + record(input: samples)`.

### Best practices for including pipelines

Pipeline inclusion only captures the pipeline's main script and included modules. It does not capture external context such as config or the `lib` directory. As a result, the pipeline should be written in a way that works when included in another pipeline:

1. Pipeline parameters should be defined in the script `params` block and referenced only in the entry workflow and `output` block. The config should only declare *config params* (params that only affect config settings).

2. Project-level assets (`projectDir`, `bin`, `lib`) should not be used since the meta-pipeline will have a different project root. Module-level assets can be safely used through the module `resources/` bundle and `moduleDir`.

3. Default `ext` settings should be specified in the process definition or avoided in favor of process inputs.

4. Software dependencies (`container`, `conda`) should be declared in the process definition, not in config.

5. Workflow outputs should be published using the `output` block, not `publishDir`.

None of these constraints are absolute. All of them can be circumvented by manually replicating the external context in the meta-pipeline. Following these constraints simply makes it easier to import a pipeline with minimal extra work.

## Open Questions

### Pipeline registry and CLI

A pipeline registry would enable remote pipeline inclusion:

```groovy
// module
include { BWA_MEM } from 'nf-core/bwa/mem'

// pipeline
include { workflow as NFCORE_RNASEQ } from 'nf-core/rnaseq'
```

This would require a `nextflow pipeline` command group for publishing and installing pipelines, similar to modules.

This effort is deferred, since it is not required for pipeline composition. Users can already compose pipelines by cloning or git submodules.

### Using plugin functions in included pipeline

If an included pipeline uses plugins, these plugins must be explicitly declared in the meta-pipeline config since they cannot be inferred from the pipeline inclusion.

If we introduce a pipeline spec, these plugin dependencies could be specified there under `requires.plugins`. When installing a pipeline, Nextflow could copy these plugin declarations into the meta-pipeline config and/or spec.

## Alternatives

### Pipeline chaining

An alternative to pipeline composition is a *pipeline chain*, in which multiple Nextflow pipelines are called in sequence via `nextflow run`.

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

- It forfeits dataflow composition. The developer must serialize/deserialize samplesheet files instead of passing channels directly between pipelines. Each pipeline must complete before the next pipeline can start.

- It requires an external workflow system instead of reusing the language that pipeline developers already know. Even the Nextflow-in-Nextflow approach shown above requires many tricks to orchestrate nested pipeline runs via the `NEXTFLOW_RUN` process.

Pipeline chaining can be practical for certain use cases, such as simple chains (A -> B -> C) of off-the-shelf pipelines. But the general solution is to compose pipelines using dataflow logic, just like any other Nextflow pipeline.

## Links

- Community issues: [#6474](https://github.com/nextflow-io/nextflow/issues/6474)
- Related: [Workflow params](20250825-workflow-params.md)
- Related: [Workflow outputs](20251020-workflow-outputs.md)
- Related: [Module system](20251114-module-system.md)
- Related: [Workflow modules](20260608-workflow-modules.md)

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
        │   ├── nextflow.config
        │   └── modules/
        └── rnaseq/
            ├── main.nf
            ├── nextflow.config
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
    input: Channel<Sample> // samplesheet
    aligner: String = 'star_salmon'
    fasta: Path
}
workflow {
    main:
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

The meta-pipeline includes each pipeline, along with its `params` block as a record type, and composes them into a new entry workflow:

```groovy
include {
    params as FetchngsParams ;
    workflow as NFCORE_FETCHNGS
} from './pipelines/nf-core/fetchngs'

include {
    params as RnaseqParams ;
    workflow as NFCORE_RNASEQ
} from './pipelines/nf-core/rnaseq'

params {
    fetchngs: FetchngsParams        // input
    strandedness: String = 'auto'   // unique to meta-pipeline
    rnaseq: RnaseqParams            // input, aligner, fasta
}

workflow {
    main:
    // fetch FASTQ samples from NCBI SRA
    samples = NFCORE_FETCHNGS( params.fetchngs )

    // adapt fetchngs output to rnaseq input (add strandedness)
    ch_samples = samples.map { r ->
        r + record(strandedness: params.strandedness)
    }

    // perform RNAseq analysis (ch_samples overrides params.rnaseq.input)
    rnaseq = NFCORE_RNASEQ( params.rnaseq + record(input: ch_samples) )

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

- **The handoff is a channel, not a file.** rnaseq declares its samplesheet input as `Channel<Sample>` instead of `Path`, so that it can be executed directly from a CSV samplesheet or called by a meta-pipeline with a live channel. When rnaseq is launched directly, the `Channel<Sample>` param is loaded from the samplesheet given on the command line. It allows rnaseq to begin aligning each sample as soon as it is emitted by fetchngs, whereas a pipeline chain would block until fetchngs finished completely.

- **Included params are partial record types.** All fields of `FetchngsParams` and `RnaseqParams` are nullable, and defaults are not pre-filled. A param that is not set (e.g. `params.rnaseq.aligner`) is `null` in the meta-pipeline, and the included pipeline applies its own default when it is called. The user can provide any rnaseq param as `--rnaseq.<name>`, the meta-pipeline can override specific params (`params.rnaseq + record(input: ch_samples)`), and the `NFCORE_RNASEQ()` call validates that all required params are present. This way, the developer only needs to declare one param for each included pipeline.

- `rnaseq.input` is supplied by the dataflow, which overrides any value given by the user.

- `rnaseq.fasta` must still be provided by the user, but the error surfaces at the `NFCORE_RNASEQ()` call rather than at launch.

- **Params and outputs are not inherited.** The meta-pipeline declares its own `params` and `output` blocks and passes params explicitly to each included pipeline. The included pipelines do not contribute any of their own params or outputs. A meta-pipeline decides for itself which outputs to publish and where to publish them.

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

In practice, the meta-pipeline will likely need to recreate the configuration shell used by the inner pipelines:

- Config params (`outdir`, `publish_dir_mode`, `max_cpus`, etc)
- Resource settings (`cpus`, `memory`, `time`, etc)
- Environment profiles (executors, software dependencies, test profiles)
- Reports (execution, timeline, trace)
- Manifest (name, authors, description, etc)
- Plugins

The meta-pipeline can reuse config files by including them, which is useful for process config. Process selectors should be written in a way that is correct both when a pipeline is executed directly *and* when it is called by a meta-pipeline.

**Remote inclusion**

With a pipeline registry (see [Pipeline registry and CLI](#pipeline-registry-and-cli)), the included pipelines would no longer need to be vendored under `pipelines/`. The includes would refer to the remote pipeline instead of the local path:

```groovy
include {
    params as FetchngsParams ;
    workflow as NFCORE_FETCHNGS
} from 'nf-core/fetchngs'

include {
    params as RnaseqParams ;
    workflow as NFCORE_RNASEQ
} from 'nf-core/rnaseq'
```

Otherwise, the meta-pipeline would work the same way.
