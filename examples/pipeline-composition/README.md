# Pipeline composition

A meta-pipeline that fetches FASTQ samples with `nf-core/fetchngs` and analyzes them with `nf-core/rnaseq`, composing the two pipelines with regular dataflow logic.

The tools are fake. Every process writes a dummy file, so the example runs anywhere in a couple of seconds.

## What it demonstrates

A pipeline (a `params` block, an entry workflow, and an `output` block) can be included as a named workflow and called like any other workflow:

```nextflow
include { workflow as NFCORE_FETCHNGS } from './pipelines/nf-core/fetchngs'

workflow {
    main:
    fetchngs = NFCORE_FETCHNGS( ids: params.fetchngs )
}
```

The `params` block acts as the `take:` section, so the pipeline is called with named arguments and params with a default value can be omitted. The `output` block acts as the `emit:` section, so `fetchngs.samples` is a channel that the calling workflow can operate on.

The two blocks can also be imported as record types, which is what `main.nf` does. That way the meta-pipeline declares one param per included pipeline instead of replicating every param:

```nextflow
include {
    params as RnaseqParams ;
    workflow as NFCORE_RNASEQ ;
    output as RnaseqOutput
} from './pipelines/nf-core/rnaseq'

params {
    rnaseq: RnaseqParams
}
```

## Layout

```
pipeline-composition/
├── main.nf                 # the meta-pipeline
├── nextflow.config
├── data/                   # accession list and reference genome
└── pipelines/
    └── nf-core/
        ├── fetchngs/main.nf
        └── rnaseq/main.nf
```

Each included pipeline is an ordinary Nextflow pipeline, vendored into the meta-pipeline repository. It can still be run on its own:

```console
$ nextflow run ./pipelines/nf-core/fetchngs --ids ./data/ids.txt
```

## Running it

```console
$ nextflow run .
```

## What to look at

**The handoff is a channel, not a file.** `rnaseq` declares its samplesheet input as `Channel<Sample>`, so `main.nf` passes it a live channel instead of a CSV file:

```nextflow
ch_samples = fetchngs.samples.map { sample ->
    sample + record(strandedness: params.strandedness)
}

rnaseq = NFCORE_RNASEQ( params.rnaseq + record(input: ch_samples) )
```

Each sample starts aligning as soon as fetchngs emits it. A pipeline chain that shells out to `nextflow run` would wait for fetchngs to finish first.

**Params flow through the imported record.** Every field of `RnaseqParams` is nullable, so the user fills in what they want:

```console
$ nextflow run . --rnaseq.aligner hisat2
```

`nextflow.config` sets `params.rnaseq.fasta`, and the two are merged. `aligner` is left unset by both, so rnaseq applies its own default. `input` is overridden by the dataflow, whatever the user passes.

**Outputs flow through the imported record.** `output { rnaseq: RnaseqOutput {} }` re-declares each rnaseq output with the directives rnaseq gave it, which is why the BAM files land in `bams/` without `main.nf` saying so. It is all-or-nothing: to publish only some of them, declare each output separately.

**Only the calling pipeline publishes.** `fetchngs` publishes its samples to `fastq/`, but nothing lands there when it is included, because `main.nf` doesn't declare that output. The outputs of an included pipeline are emitted to the calling workflow, which decides what to publish.

**Processes are scoped by the include alias.** The run log shows `NFCORE_RNASEQ:STAR_ALIGN`, and `nextflow.config` targets it the same way:

```nextflow
process {
    withName: 'NFCORE_RNASEQ:STAR_ALIGN' {
        cpus   = 2
        memory = 2.GB
    }
}
```

## See also

- [Pipeline composition](../../adr/20260608-pipeline-composition.md): the design decision behind this example.
- [Typed workflows](../../docs/workflow-typed.mdx): the syntax used here.
