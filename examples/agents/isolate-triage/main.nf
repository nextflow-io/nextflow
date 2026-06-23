nextflow.enable.types = true

/*
 * ============================================================================
 *  REAL-WORLD USE CASE: autonomous bacterial-isolate triage
 * ============================================================================
 *
 *  A surveillance / clinical-micro lab receives short reads from a bacterial
 *  isolate. The agent is asked, in plain language, to assemble the genome and
 *  decide whether the assembly is good enough to annotate — and only then run
 *  the (expensive) annotation step. This is a DATA-DRIVEN, adaptive pipeline:
 *  the path through the tools is decided at run time from intermediate results,
 *  not hard-wired in the DAG.
 *
 *  Three nf-core modules are included so they are discoverable by `module_run`:
 *    1. nf-core/skesa        — assemble reads -> contigs (a FASTA file handle)
 *    2. nf-core/assemblyscan — compute assembly statistics -> a small JSON file
 *    3. nf-core/prokka       — annotate contigs -> annotation files
 *
 *  TWO CAPABILITY TOOLS:
 *  -----------------------------------------------------------------------
 *  `tools 'module_run', 'filesystem'`
 *
 *  `module_run` exposes a single generic tool whose `module` enum lists every
 *  in-scope process (all three included modules above). The LLM picks which
 *  module to run and in what order — the adaptive QC gate is a reasoning step,
 *  not a hard-coded edge in the DAG.
 *
 *  `filesystem` gives the agent a sandboxed read/write/list/exists tool scoped
 *  to its per-invocation work directory (plus module-output paths). The agent
 *  uses it to write a short triage summary JSON to the work dir alongside the
 *  verdict it emits on the output channel.
 *
 *  TWO KINDS OF TOOL OUTPUT (why this works with the LLM-tool model):
 *  -------------------------------------------------------------------
 *  Bulk/binary artifacts the LLM only forwards between tools — the skesa
 *  contigs, the prokka annotations — come back as OPAQUE ABSOLUTE-PATH HANDLES:
 *  the model never reads their bytes, it just passes the path to the next tool.
 *
 *  But the QC GATE needs the model to REASON OVER numbers (N50, #contigs), and
 *  assemblyscan emits those as a SMALL `.json` file. Nextflow detects that this
 *  output is a small, text/structured file (extension in a known set, size under
 *  the `agent.maxToolOutputInlineSize` cap, not binary) and INLINES ITS CONTENTS
 *  into the tool result the LLM receives — so the model reads the real stats and
 *  gates on them. A large or binary output would stay an opaque handle. You do
 *  not annotate anything for this: it is inferred from the output file's format
 *  and size at run time.
 *
 *  SELF-DESCRIBING TOOLS: each nf-core module describes itself to the LLM from
 *  its meta.yml — the tool description and input schema (including the `meta` map
 *  and its `id`) are wired into the module_run tool spec automatically. That is
 *  why the `instruction` below is a purely HIGH-LEVEL functional goal: it never
 *  names a tool, an argument shape, or `meta.id` — the LLM learns all of that
 *  from the module schemas surfaced via `module_run`.
 */

// Include the three nf-core modules; they are discovered by `module_run` at
// agent-run time and surfaced as the `module` enum in the single tool.
include { SKESA }        from 'nf-core/skesa'
include { ASSEMBLYSCAN } from 'nf-core/assemblyscan'
include { PROKKA }       from 'nf-core/prokka'

record Isolate {
    sample_id: String
    organism:  String      // best-guess species, e.g. "Escherichia coli"
    reads:     String      // path to the short-read FASTQ (a handle for the tools)
}

agent triage {
    model 'openai/gpt-5-mini'
    // `goal` states the high-level objective; it is folded into the system message
    // alongside `instruction` to steer the multi-turn module_run loop. It is
    // advisory — `maxIterations` remains the hard cap on the tool-calling loop.
    goal 'Assemble the isolate, QC-gate the assembly, annotate only if it passes, and report a clear PASS/FAIL verdict with a short written summary.'
    // NOTE: the QC thresholds below are deliberately LENIENT so the bundled demo
    // reads (a small ggal fragment that assembles to N50≈863, ~24 contigs) PASS the
    // gate and the annotation branch (prokka) runs. For a real bacterial isolate use
    // realistic values, e.g. "N50 below 20000 bp OR more than 500 contigs".
    instruction '''\
        You triage bacterial isolate assemblies. Work step by step:
          1. Assemble the isolate's sequencing reads into contigs using module_run.
          2. Compute assembly QC statistics from the contigs using module_run.
          3. QC GATE: read the assembly statistics and decide. If the assembly is
             too fragmented — N50 below 500 bp OR more than 1000 contigs — it
             FAILS QC: write a brief summary JSON to the work directory using the
             filesystem tool, then reply
             "FAIL <sample>: fragmented (N50=<n50>, contigs=<contigs>), needs manual review"
             and DO NOT annotate.
          4. Otherwise the assembly PASSES QC: annotate the contigs using module_run,
             write a brief summary JSON to the work directory using the filesystem tool,
             then reply
             "PASS <sample>: N50=<n50>, contigs=<contigs>, annotation=<path to the annotation file>".
        '''.stripIndent()

    // `module_run` discovers all included modules (skesa, assemblyscan, prokka)
    // and surfaces them as a single tool with a `module` enum.
    // `filesystem` enables sandboxed file read/write/list/exists in the agent work dir.
    tools 'module_run', 'filesystem'

    input:
        isolate: Isolate
    output:
        verdict: String          // tools => plain (non-record) output; the LLM's final report

    prompt:
    """
    Triage isolate '${isolate.sample_id}' (${isolate.organism}).
    Short reads: ${isolate.reads}
    """
}

workflow {
    // NOTE: this example needs a real FASTQ at `data/sample.fastq` (the `data/`
    // dir is gitignored — provide your own reads). The nf-core modules run in
    // containers, so a container runtime (Docker/Wave) is also required.
    triage(channel.of(
        record(sample_id: 'isolate_001', organism: 'Escherichia coli',
               reads: "${projectDir}/data/sample.fastq")
    ))
    .view { v -> "TRIAGE: ${v}" }
}
