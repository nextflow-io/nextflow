nextflow.enable.types = true

// Real-world adaptive agent: assemble (SKESA) → QC stats (ASSEMBLYSCAN) → a
// data-driven QC gate → annotate (PROKKA) only if it passes, writing a summary via the
// `filesystem` tool. The route through the tools is decided at run time, and
// ASSEMBLYSCAN's small JSON stats are inlined into the tool result so the model can gate
// on N50/#contigs. See README.md.
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
    // `goal` = high-level objective, folded into the system message (advisory).
    goal 'Assemble the isolate, QC-gate the assembly, annotate only if it passes, and report a clear PASS/FAIL verdict with a short written summary.'
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
