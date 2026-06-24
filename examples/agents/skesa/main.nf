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
 * processes (both included modules and locally-defined ones) and surfaces EACH
 * of them to the LLM as its OWN tool, named after the module (here `SKESA`),
 * whose parameters schema IS the module's flattened input schema.
 * No tool-by-name `tools 'nf-core/skesa'` string is needed.
 */
include { SKESA } from 'nf-core/skesa'

/*
 * An agent that assembles a microbial genome by calling the nf-core/skesa
 * module as a tool via the `module_run` capability.
 *
 *   - `tools 'module_run'` exposes every in-scope process as its OWN tool. The
 *     runnable modules are discovered from the script's `include` statements
 *     and any locally-defined processes. The included process is named `SKESA`
 *     (nf-core process names are upper-case), so the LLM sees a tool named
 *     `SKESA` whose parameters schema is {meta, fastq} (required, derived from
 *     the registry/meta.yml) — the model is schema-constrained to those names.
 *
 *   - The LLM reads this agent's prompt (built from `AssemblyRequest`) and
 *     SYNTHESISES the tool call, e.g.:
 *       SKESA({"meta":{"id":"sample1"},"fastq":"/data/sample1.fastq.gz"})
 *     It maps `sample_id` -> meta.id and `reads` -> fastq; the two
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
    instruction 'Use the SKESA tool with {meta:{id:<sample_id>}, fastq:<reads path>} to assemble the provided reads, then report the path to the assembled contigs.'

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
