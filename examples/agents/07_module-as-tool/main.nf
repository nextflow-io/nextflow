nextflow.enable.types = true

// Include nf-core/skesa so `nf:module_run` surfaces it as the `SKESA` tool.
include { SKESA } from 'nf-core/skesa'

record Assembly {
    sample_id: String
    contigs: Path
}

// Agent that calls the SKESA tool to assemble reads into contigs. The input record
// is destructured, and its shape need NOT match skesa's tool input ({meta, fastq}) --
// the LLM bridges the two from the prompt. The `contigs` field of the typed output is the
// path of the contigs returned by the tool call. See README.md.
agent assembler {
    model 'openai/gpt-5-mini'
    instruction """
        You are a genome assembly assistant. Use the available tools to
        assemble the provided sequencing reads into contigs,
        then report the path to the assembled contigs.
        """

    tools 'nf:module_run'

    input:
    record(sample_id: String, reads: Path)
    output:
    assembly: Assembly

    prompt:
    """
    Assemble the genome for sample '${sample_id}'.
    The input FASTQ reads are at: ${reads}
    """
}

workflow {
    // NOTE: this example needs a FASTQ at `data/sample.fastq`. `data` is a symlink to
    // `examples/data/`, which is shared by every example needing input and fetched once --
    // see examples/data/README.md for the command.
    // nf-core/skesa runs in a container, so Docker + Wave (or another container
    // runtime) is also required.
    assembler(channel.of(
        record(sample_id: 'sample1', reads: file("${projectDir}/data/sample.fastq"))
    ))
    .view { a -> "ASSEMBLY ${a.sample_id}=${a.contigs}" }
}
