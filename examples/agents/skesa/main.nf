nextflow.enable.types = true

// The agent's input record. Its shape need NOT match skesa's tool input
// ({meta, fastq}) — the LLM bridges the two from the prompt. See README.md.
record AssemblyRequest {
    sample_id: String
    reads: String          // absolute path to a FASTQ file (an opaque path handle)
}

// Include nf-core/skesa so `module_run` surfaces it as the `SKESA` tool.
include { SKESA } from 'nf-core/skesa'

// Agent that calls the SKESA tool to assemble reads into contigs. An agent that
// declares `tools` must use a plain output type (here `Path`), not a record.
agent assembler {
    model 'openai/gpt-5-mini'
    instruction 'Use the SKESA tool to assemble the provided reads, then report the path to the assembled contigs.'

    tools 'module_run'

    input:
        req: AssemblyRequest
    output:
        assembly_path: Path

    prompt:
    """
    Assemble the genome for sample '${req.sample_id}'.
    The input FASTQ reads are at: ${req.reads}
    """
}

workflow {
    // NOTE: this example needs a FASTQ at `data/sample.fastq` (the `data/` dir is
    // gitignored). Fetch the sarscov2 test dataset used for demos:
    //   mkdir -p data && curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz | gunzip > data/sample.fastq
    // nf-core/skesa runs in a container, so Docker + Wave (or another container
    // runtime) is also required.
    assembler(channel.of(
        record(sample_id: 'sample1', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { path -> "ASSEMBLY=${path}" }
}
