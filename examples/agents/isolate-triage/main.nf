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
 *  `module_run` exposes EACH in-scope process as its OWN tool (all three
 *  included modules above become the tools SKESA, ASSEMBLYSCAN, PROKKA), each
 *  with its own enforced input schema. The LLM picks which module to run and in
 *  what order — the adaptive QC gate is a reasoning step, not a hard-coded edge
 *  in the DAG.
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
 *  its registry/meta.yml metadata — the tool description and input schema
 *  (including the `meta` map and its `id`) become that module's OWN tool's
 *  enforced parameters schema automatically. That is why the `instruction`
 *  below is a purely HIGH-LEVEL functional goal: it never names an argument
 *  shape or `meta.id` — the LLM learns all of that from the per-module tool
 *  schemas surfaced via `module_run`.
 */

// Include the three nf-core modules; they are discovered by `module_run` at
// agent-run time and each surfaced as its OWN tool (SKESA, ASSEMBLYSCAN, PROKKA).
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
    // NOTE: the QC thresholds below are deliberately LENIENT (N50 < 500 bp OR
    // more than 1000 contigs). With the sarscov2 test FASTQ (see workflow block
    // below) the assembly reaches ~N50=310, ~6 contigs and CORRECTLY FAILS this
    // gate — the agent returns a FAIL verdict and skips PROKKA annotation.
    // To exercise the PASS + PROKKA branch, use a real bacterial isolate FASTQ
    // (better-assembling reads) or tighten the thresholds to match the test data.
    instruction '''\
        You triage bacterial isolate assemblies. Work step by step:
          1. Assemble the isolate's sequencing reads into contigs using the SKESA tool.
          2. Compute assembly QC statistics from the contigs using the ASSEMBLYSCAN tool.
          3. QC GATE: read the assembly statistics and decide. If the assembly is
             too fragmented — N50 below 500 bp OR more than 1000 contigs — it
             FAILS QC: write a brief summary JSON to the work directory using the
             filesystem tool, then reply
             "FAIL <sample>: fragmented (N50=<n50>, contigs=<contigs>), needs manual review"
             and DO NOT annotate.
          4. Otherwise the assembly PASSES QC: annotate the contigs using the PROKKA tool,
             write a brief summary JSON to the work directory using the filesystem tool,
             then reply
             "PASS <sample>: N50=<n50>, contigs=<contigs>, annotation=<path to the annotation file>".
        '''.stripIndent()

    // `module_run` discovers all included modules (skesa, assemblyscan, prokka)
    // and surfaces EACH as its OWN tool (SKESA, ASSEMBLYSCAN, PROKKA) with an
    // enforced input schema.
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
    // NOTE: this example needs a FASTQ at `data/sample.fastq` (the `data/` dir is
    // gitignored). Fetch the sarscov2 test dataset used here:
    //   mkdir -p data && curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz | gunzip > data/sample.fastq
    // The nf-core modules run in containers, so Docker + Wave (or another container
    // runtime) is also required.
    //
    // EXPECTED BEHAVIOUR with the sarscov2 test FASTQ: the assembly assembles to
    // ~N50=310, ~6 contigs — which FAILS the lenient gate (N50 < 500), so the agent
    // returns a FAIL verdict and skips PROKKA annotation. That is the correct and
    // intended gate behaviour for this small viral read set. To exercise the PASS +
    // PROKKA branch, use a real bacterial isolate FASTQ or adjust the N50/contig
    // thresholds in the `instruction` above.
    triage(channel.of(
        record(sample_id: 'isolate_001', organism: 'Escherichia coli',
               reads: "${projectDir}/data/sample.fastq")
    ))
    .view { v -> "TRIAGE: ${v}" }
}
