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
 * Include the nf-core/skesa module so its process is in scope. Under the
 * capability model the `tools 'module_run'` directive discovers all in-scope
 * processes (both included modules and locally-defined ones) and surfaces them
 * to the LLM as a single `module_run(module, args)` tool with a `module` enum.
 * No tool-by-name `tools 'nf-core/skesa'` string is needed.
 */
include { skesa } from 'nf-core/skesa'

/*
 * An agent that assembles a microbial genome by calling the nf-core/skesa
 * module as a tool via the `module_run` capability.
 *
 *   - `tools 'module_run'` enables the generic `module_run(module, args)` tool.
 *     The runnable modules are discovered from the script's `include` statements
 *     and any locally-defined processes. Here `skesa` is the sole included module,
 *     so the LLM's `module` enum has one entry.
 *
 *   - The LLM reads this agent's prompt (built from `AssemblyRequest`) and
 *     SYNTHESISES the module_run call, e.g.:
 *       module_run({"module":"skesa","args":{"meta":{"id":"sample1"},"fastq":"/data/sample1.fastq.gz"}})
 *     It maps `sample_id` -> args.meta.id and `reads` -> args.fastq; the two
 *     structures never have to match.
 *
 *   - skesa runs as a real dataflow node (its container, work dir, caching), and
 *     its `fasta` output { meta, *.fa } comes back to the LLM as JSON with the
 *     FASTA path as an absolute path handle.
 *
 * NOTE: an agent that declares `tools` must use a PLAIN output type (here `Path`)
 * — combining tools with a structured record output is not supported in v1.
 */
agent assembler {
    model 'openai/gpt-5-mini'
    instruction 'Use the module_run tool to run the skesa module on the provided reads, then report the path to the assembled contigs.'

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
    assembler(channel.of(
        record(sample_id: 'sample1', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { path -> "ASSEMBLY=${path}" }
}
