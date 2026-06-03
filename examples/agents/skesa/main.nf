nextflow.enable.types = true

/*
 * The AGENT's input record. Note its shape is NOTHING like skesa's tool input
 * (which is a tuple of { meta: map, fastq: file }). There is no field-to-field
 * mapping between this record and the tool — the LLM bridges the two.
 */
record AssemblyRequest {
    sample_id: String
    reads: String          // absolute path to a FASTQ file (an opaque path handle)
}

/*
 * An agent that assembles a microbial genome by calling the real nf-core/skesa
 * module as a tool.
 *
 *   - `tools 'nf-core/skesa'` resolves the module from the registry. Its tool
 *     input schema is derived from skesa's meta.yml: the input tuple is FLATTENED
 *     to top-level properties { meta (object), fastq (string/path) }, each carrying
 *     skesa's own field descriptions — the LLM gets a documented tool for free.
 *
 *   - The LLM reads this agent's prompt (built from `AssemblyRequest`) and SYNTHESIZES
 *     the skesa call, e.g. skesa({ "meta": {"id": "sample1"}, "fastq": "/data/sample1.fastq.gz" }).
 *     It maps `sample_id` -> meta.id and `reads` -> fastq itself; the two structures
 *     never have to match.
 *
 *   - skesa runs as a real dataflow node (its container, work dir, caching), and its
 *     `fasta` output { meta, *.fa } comes back to the LLM as JSON with the FASTA path
 *     as an absolute path handle.
 *
 * NOTE: an agent that declares `tools` must use a PLAIN output type (here `String`) —
 * combining tools with a structured record output is not supported in v1.
 */
agent assembler {
    model 'openai/gpt-5-mini'
    instruction 'Assemble the genome from the provided sequencing reads and report the path to the assembled contigs.'

    tools 'nf-core/skesa'

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
    assembler(channel.of(
        record(sample_id: 'sample1', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { path -> "ASSEMBLY=${path}" }
}
