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
 *  The agent has three tools:
 *    1. nf-core/skesa   — assemble reads -> contigs (a FILE handle)
 *    2. assembly_qc     — compute assembly metrics -> SCALAR values
 *    3. nf-core/prokka  — annotate contigs -> annotation FILES
 *
 *  KEY DESIGN PRINCIPLE (why this works with the LLM-tool model):
 *  ----------------------------------------------------------------
 *  The LLM never reads file CONTENTS — it only sees opaque absolute-path
 *  handles and passes them between tools. So anything the LLM must REASON OVER
 *  (here: "is this assembly good enough?") must be returned as a SCALAR / `val`
 *  tool output, not as a file. That is exactly what `assembly_qc` does: it turns
 *  a FASTA file into numbers (N50, #contigs) the model can compare to a
 *  threshold. Heavy artifacts (contigs, GFF/GBK) stay as handles the model just
 *  forwards. Design your tool outputs around this split.
 */

record Isolate {
    sample_id: String
    organism:  String      // best-guess species, e.g. "Escherichia coli"
    reads:     String      // path to the short-read FASTQ (a handle for the tools)
}

/*
 * A thin, CONTAINER-FREE QC tool: it parses the assembled FASTA and emits
 * decision-relevant metrics as SCALARS the LLM can reason over. Pure Groovy
 * `exec:` — no external binary, so the model gets numbers, not a report file.
 */
process assembly_qc {
    input:
        contigs: Path
    output:
        n50:         Integer
        num_contigs: Integer
        total_bp:    Integer
    exec:
        final lengths = []
        int cur = 0
        contigs.eachLine { String line ->
            if( line.startsWith('>') ) { if( cur > 0 ) lengths.add(cur); cur = 0 }
            else cur += line.trim().length()
        }
        if( cur > 0 ) lengths.add(cur)
        lengths.sort()
        num_contigs = lengths.size()
        total_bp = (lengths.sum() ?: 0) as Integer
        // N50: the contig length at which half the assembly is in contigs >= that length
        final half = total_bp.intdiv(2)
        int acc = 0; int n = 0
        for( int i = lengths.size() - 1; i >= 0; i-- ) {
            acc += lengths[i]
            if( acc >= half ) { n = lengths[i]; break }
        }
        n50 = n
}

agent triage {
    model 'openai/gpt-5-mini'
    instruction '''\
        You triage bacterial isolate assemblies. Work step by step:
          1. Assemble the isolate's sequencing reads into contigs.
          2. Compute assembly QC metrics (n50, num_contigs, total_bp) from the contigs.
          3. QC GATE: if n50 < 20000 OR num_contigs > 500, the assembly is too fragmented —
             reply "FAIL <sample>: fragmented (N50=<n50>, contigs=<num_contigs>), needs manual review"
             and DO NOT annotate.
          4. Otherwise, annotate the contigs, then reply
             "PASS <sample>: N50=<n50>, contigs=<num_contigs>, annotation=<path to the .gff>".
        '''.stripIndent()

    // `tools` is varargs and freely mixes registry refs with in-scope process names:
    tools 'nf-core/skesa', 'assembly_qc', 'nf-core/prokka'

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
    triage(Channel.of(
        record(sample_id: 'isolate_001', organism: 'Escherichia coli',
               reads: "${projectDir}/data/sample.fastq")
    ))
    .view { v -> "TRIAGE: ${v}" }
}
