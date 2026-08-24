nextflow.enable.types = true

// Runner-native tools (`fs:*`, `shell:bash`), executed by the runner rather than brokered back to
// the driver. `shell:bash` is `pi`-only. See README.md.

// Every field is checkable against the file, which makes this example a test, not a demo.
record AssemblyStats {
    contig_count: Integer
    total_bases: Integer
    gc_percent: Float
    n50: Integer
    longest_contig: Integer
}

agent assembly_qc {
    model 'openai/gpt-5-mini'

    instruction '''\
        You are an assembly QC assistant with a POSIX shell in your working directory.

        The FASTA is already staged in your working directory under the name the prompt gives you.
        Compute every number from that file with the shell — `awk`, `grep`, `sort`, `wc`.

        Do NOT count by reading. Reading thousands of bases and estimating is simply a wrong
        answer.

        The FASTA is line-wrapped: a contig's sequence is EVERY line after its `>` header up to the
        next header, so concatenate those lines before measuring one. A wrapped line is not a
        contig — if your longest contig comes out equal to the wrap width, you measured lines.

        Definitions, so the arithmetic is unambiguous:
          - contig_count    number of `>` header lines
          - total_bases     total sequence characters, headers and newlines excluded
          - gc_percent      100 * (G + C) / total_bases, case-insensitive, 2 decimal places
          - n50             per-contig lengths sorted descending, accumulated; report the length of
                            the contig at which the running total first reaches half of total_bases
          - longest_contig  length of the longest contig, concatenated across its wrapped lines

        Write your working to `stats.tsv` so the numbers can be audited, then report them.
        '''.stripIndent()

    tools 'fs:*', 'shell:bash'

    input:
    contigs: Path

    output:
    stats: AssemblyStats

    prompt:
    """
    Compute the assembly statistics for the FASTA file ${contigs} in your working directory.
    """
}

workflow {
    // A `Path` input is staged into the agent's task directory and bind-mounted into the runner
    // container, so the agent's shell opens it under its plain name — see README.md.
    assembly_qc(channel.fromPath("${moduleDir}/data/contigs.fa"))
        .view { s ->
            """\
            contigs   = ${s.contig_count}
            bases     = ${s.total_bases}
            GC%       = ${s.gc_percent}
            N50       = ${s.n50}
            longest   = ${s.longest_contig}
            """.stripIndent()
        }
}
