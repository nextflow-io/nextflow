nextflow.enable.types = true

// The agent's input record.
record Sample {
    sample_id: String
    reads: String          // absolute path to a FASTQ file (an opaque path handle)
}

// Include the nf-core modules; `module_run` surfaces each as a tool (SKESA, ASSEMBLYSCAN).
include { SKESA }        from 'nf-core/skesa'
include { ASSEMBLYSCAN } from 'nf-core/assemblyscan'

// Goal-directed tool composition: from a high-level `goal` (no step list) the model
// plans assemble (SKESA) → assembly-stats (ASSEMBLYSCAN) → judge against an N50 bar.
// This is composition, not iteration — contrast with `convergence-loop`. See README.md.
agent qc {
    model 'openai/gpt-5-mini'

    // `instruction` = the role + constraints (HOW the agent may act). It does
    // NOT list the steps and does NOT name any tool's input fields — those come
    // from the goal and the registry-derived tool schemas.
    instruction '''\
        You are a genome-assembly QC assistant. Use the available tools to do
        the work; pick whichever tools fit each step. When a tool returns a path
        to a results file, read that file to obtain the actual numbers. Base
        every metric on tool output — never guess a value.
        '''.stripIndent()

    // `goal` = the OBJECTIVE that drives the multi-turn loop. The model keeps
    // calling tools until it has what the goal asks for, then stops.
    goal '''\
        Assemble the provided sequencing reads into contigs, compute the
        assembly quality statistics (N50 and the number of contigs) from those
        contigs, then report whether the assembly meets the quality bar of
        N50 >= 500 bp. Include the N50, the contig count, and the absolute path
        to the contigs FASTA in your final answer.
        '''.stripIndent()

    tools 'module_run', 'filesystem'
    maxIterations 15

    input:
        sample: Sample
    output:
        report: String       // tools => plain (non-record) output; the LLM's final report

    prompt:
    """
    Assess the assembly quality for sample '${sample.sample_id}'.
    The sequencing reads (FASTQ) are at: ${sample.reads}
    """
}

workflow {
    // NOTE: this example needs a FASTQ at `data/sample.fastq` (the `data/` dir is
    // gitignored). Fetch the sarscov2 test dataset used for demos:
    //   mkdir -p data && curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz | gunzip > data/sample.fastq
    // The nf-core modules run in containers, so Docker + Wave (or another
    // container runtime) is also required.
    //
    // EXPECTED with the sarscov2 test FASTQ: it assembles to ~N50=310, ~6
    // contigs, which is BELOW the N50>=500 bar, so the agent reports the
    // assembly does NOT meet the quality bar. That is the correct, deterministic
    // outcome for this small viral read set — the point of the demo is the
    // goal-directed loop, not the verdict.
    qc(channel.of(
        record(sample_id: 'sample1', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { r -> "QC=${r}" }
}
