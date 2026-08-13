nextflow.enable.types = true

// A tool agent fanning out over a queue in PARALLEL, paired with structured output.
// See README.md.

// One sequence's annotation, produced by the agent (structured output).
record SeqAnnotation {
    seq:     String   // the original sequence (upper-cased, validated)
    revcomp: String   // the reverse complement — MUST come from the `revcomp` tool
    length:  Integer  // sequence length
    note:    String   // any caveat (e.g. non-ACGT bases seen) or "ok"
}

// Deterministic tool: reverse-complement a DNA string. `nf:module_run` exposes this in-scope
// process to the LLM as the tool `revcomp` (input seq:String -> output rc:String).
process revcomp {
    input:
    seq: String
    output:
    rc: String
    exec:
        def comp = [A:'T', T:'A', C:'G', G:'C', N:'N']
        rc = seq.toUpperCase().reverse().collect { comp[it as String] ?: 'N' }.join()
}

// The agent: for each sequence it validates the input, calls the `revcomp` tool for the
// exact reverse complement, and returns a structured SeqAnnotation. Combining `tools` with
// a record output runs the tool loop schema-free then a final structuring turn encodes the
// answer into the SeqAnnotation schema (M5), now on the task path.
agent annotate {
    model 'openai/gpt-5-mini'
    instruction '''\
        You annotate short DNA sequences. For the input sequence:
          - Upper-case it and report its length.
          - To obtain the reverse complement you MUST call the `revcomp` tool with the
            sequence; never compute the reverse complement yourself.
          - If the sequence contains any base other than A/C/G/T, say so in `note`
            (otherwise set `note` to "ok").
        '''.stripIndent()

    tools 'nf:module_run'

    input:
    seq: String
    output:
    annotation: SeqAnnotation

    prompt:
    """
    Annotate this DNA sequence:

    ${seq}
    """
}

workflow {
    // A small queue of sequences (one has an ambiguous base to exercise the `note` path).
    // Each item runs `annotate` as its own parallel agent TaskRun, each making a `revcomp`
    // tool call — the whole point of tools-on-the-task-path.
    def seqs = channel.of(
        'ATGCATTAGCCG',
        'GGGGCCCCTTTT',
        'ACGTACGTNACG',   // contains an N -> agent should flag it in `note`
        'TTACGGATCCAA',
    )

    annotate(seqs).view { a ->
        """\
        len=${a.length}  note=${a.note}
           seq     = ${a.seq}
           revcomp = ${a.revcomp}
        """.stripIndent()
    }
}
