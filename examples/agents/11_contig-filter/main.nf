nextflow.enable.types = true

// Convergence loop that drives REAL nf-core modules many times. The agent titrates
// sequencing depth: it re-runs SEQTK_SAMPLE (subsample the reads) + SKESA (assemble)
// with a varying `sample_size`, measures the resulting N50 each round, and converges
// on the SMALLEST subsample whose assembly still clears the N50 bar. `sample_size` is
// a declared module input, so `nf:module_run` exposes it as a tunable tool parameter —
// the knob the loop turns. Contrast `convergence-loop` (a local tool). See README.md.

include { SEQTK_SAMPLE } from 'nf-core/seqtk/sample'
include { SKESA }        from 'nf-core/skesa'

record Isolate {
    sample_id: String
    reads: String          // absolute path to a FASTQ file (an opaque path handle)
}

// Local tool `assembly_stats`: compute N50 and contig count from a contigs FASTA.
// (A cheap exec: step so measuring the metric doesn't spin a container each round.)
process assembly_stats {
    input:
        // Stays `String`: a `Path` tool input is rejected by the agent tool schema
        // ("not yet supported as an agent tool"), so the agent hands the path over as text.
        contigs: String
    output:
    result: String
    exec:
        def lengths = []
        def cur = 0
        // `file(...)`, not `new File(...)`: the assembler emits an `s3://`/`az://`/`gs://`
        // URI whenever the run has a remote work dir, and only Nextflow's `file()` resolves
        // those through the matching filesystem provider. `new File()` silently mangles a
        // URI into a relative local path and fails with FileNotFoundException.
        file(contigs).eachLine { line, n ->
            def s = line.trim()
            if( s.startsWith('>') ) {
                if( cur > 0 ) { lengths.add(cur) }
                cur = 0
            }
            else {
                cur = cur + s.length()
            }
        }
        if( cur > 0 ) { lengths.add(cur) }

        def sorted = lengths.sort().reverse()          // descending
        def total = 0
        sorted.each { L -> total = total + L }

        def half = total / 2.0
        def cum = 0
        def n50 = 0
        sorted.each { L ->
            cum = cum + L
            if( n50 == 0 && cum >= half ) { n50 = L }
        }

        result = String.format('contigs=%d total_length=%d n50=%d',
                               sorted.size(), total as Integer, n50 as Integer)
}

agent optimizer {
    model 'openai/gpt-5-mini'

    // role/constraints only — no tool-field structure, no fixed step list.
    instruction '''\
        You right-size sequencing depth for genome assembly. Use the available
        tools to subsample the reads to a given size, assemble the subsample, and
        measure the assembly. The tools are the ONLY way to learn the N50 for a
        given subsample — never guess, always call them. Search the subsample size
        coarsely, then refine.
        '''.stripIndent()

    // the OBJECTIVE that drives the convergence loop.
    goal '''\
        Find the SMALLEST read subsample that still assembles well. For a candidate
        subsample size (a fraction between 0 and 1 of the reads), subsample the
        reads, assemble the subsample into contigs, and measure the N50. Smaller
        subsamples generally give a lower (worse) N50, so binary-search the
        fraction: start at 0.5, then halve or raise the search toward the smallest
        fraction whose assembly still has N50 >= 300 bp. Evaluate at most 6
        fractions — do NOT scan exhaustively. Report the optimal subsample size,
        the resulting N50, the contig count, and the path to the contigs FASTA.
        '''.stripIndent()

    tools 'nf:module_run'
    maxIterations 25

    input:
    isolate: Isolate
    output:
    report: String

    prompt:
    """
    Right-size the sequencing depth for isolate '${isolate.sample_id}'.
    Sequencing reads (FASTQ): ${isolate.reads}
    """
}

workflow {
    // NOTE: needs the GZIPPED FASTQ at `data/sample.fastq.gz` -- SEQTK_SAMPLE names its
    // output after the input and emits `*.fastq.gz`, so the input must be gzipped, unlike
    // the sibling examples that read `sample.fastq`. `data` is a symlink to
    // `examples/data/`, which holds both forms and is fetched once -- see
    // examples/data/README.md for the command.
    // SEQTK_SAMPLE and SKESA run in containers, so Docker + Wave (or another
    // runtime) is required — and each loop iteration runs both, so this is slower
    // than the local-tool `convergence-loop`.
    optimizer(channel.of(
        record(sample_id: 'isolate_001', reads: "${projectDir}/data/sample.fastq.gz")
    ))
    .view { r -> "RESULT=${r}" }
}
