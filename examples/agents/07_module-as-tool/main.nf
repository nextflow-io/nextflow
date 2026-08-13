nextflow.enable.types = true

// The agent's input record. Its shape need NOT match skesa's tool input
// ({meta, fastq}) — the LLM bridges the two from the prompt. See README.md.
record AssemblyRequest {
    sample_id: String
    reads: Path          // absolute path to a FASTQ file (an opaque path handle)
}

// Include nf-core/skesa so `nf:module_run` surfaces it as the `SKESA` tool.
include { SKESA } from 'nf-core/skesa'

// Agent that calls the SKESA tool to assemble reads into contigs. An agent that
// declares `tools` must use a plain output type (here `Path`), not a record.
agent assembler {
    model 'openai/gpt-5-mini'
    instruction """
        You are a genome assembly assistant. Use the available tools to
        assemble the provided sequencing reads into contigs,
        then report the path to the assembled contigs.
        """

    tools 'nf:module_run'

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
    // NOTE: this example needs a FASTQ at `data/sample.fastq`. `data` is a symlink to
    // `examples/data/`, which is shared by every example needing input and fetched once --
    // see examples/data/README.md for the command.
    // nf-core/skesa runs in a container, so Docker + Wave (or another container
    // runtime) is also required.
    assembler(channel.of(
        record(sample_id: 'sample1', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { path -> "ASSEMBLY=${path}" }
}
