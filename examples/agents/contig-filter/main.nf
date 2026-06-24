nextflow.enable.types = true

// Convergence loop on real module output: SKESA assembles once, then the agent re-runs
// the SAME filter_contigs tool many times, varying `min_len`, to maximise N50 subject to
// a retained-length floor. `min_len` is a declared input. See README.md.

include { SKESA } from 'nf-core/skesa'

record Isolate {
    sample_id: String
    reads: String          // absolute path to a FASTQ file (an opaque path handle)
}

// Tool `filter_contigs`: drop contigs shorter than min_len; report N50/retained/%/counts.
process filter_contigs {
    input:
        contigs: String
        min_len: Integer
    output:
        result: String
    exec:
        def lengths = []
        def cur = 0
        new File(contigs).eachLine { line, n ->
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

        def totalAll = 0
        lengths.each { L -> totalAll = totalAll + L }

        def kept = []
        lengths.each { L -> if( L >= min_len ) { kept.add(L) } }
        kept = kept.sort().reverse()                  // descending

        def retained = 0
        kept.each { L -> retained = retained + L }
        def pct = totalAll > 0 ? (retained / totalAll) : 0.0

        def half = retained / 2.0
        def cum = 0
        def n50 = 0
        kept.each { L ->
            cum = cum + L
            if( n50 == 0 && cum >= half ) { n50 = L }
        }

        result = String.format('min_len=%d kept_contigs=%d total_contigs=%d retained_length=%d total_length=%d pct_retained=%.3f n50=%d',
                               min_len as Integer, kept.size(), lengths.size(), retained as Integer, totalAll as Integer, pct as Double, n50 as Integer)
}

agent optimizer {
    model 'openai/gpt-5-mini'

    // role/constraints only — no tool-field structure, no fixed step list.
    instruction '''\
        You improve a genome assembly by tuning a contig length filter. First
        assemble the reads to get contigs, then use the filter_contigs tool to
        measure the effect of a minimum contig length. filter_contigs is the ONLY
        way to learn the N50 and retained fraction for a given min_len — never
        guess, always call it. Search min_len coarsely, then refine.
        '''.stripIndent()

    // the OBJECTIVE that drives the convergence loop.
    goal '''\
        Produce the best-quality assembly you can: assemble the reads, then find
        the minimum contig length that MAXIMISES N50 while RETAINING at least 70%
        of the total assembly length. Iterate min_len until N50 stops improving
        within that constraint. Report the optimal min_len, the resulting N50, the
        percent retained, and the path to the contigs.
        '''.stripIndent()

    tools 'module_run'
    maxIterations 30

    input:
        isolate: Isolate
    output:
        report: String

    prompt:
    """
    Optimise the assembly for isolate '${isolate.sample_id}'.
    Sequencing reads (FASTQ): ${isolate.reads}
    """
}

workflow {
    // NOTE: needs a FASTQ at `data/sample.fastq` (the `data/` dir is gitignored).
    // Fetch the sarscov2 test dataset:
    //   mkdir -p data && curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz | gunzip > data/sample.fastq
    // SKESA runs in a container, so Docker + Wave (or another runtime) is required.
    optimizer(channel.of(
        record(sample_id: 'isolate_001', reads: "${projectDir}/data/sample.fastq")
    ))
    .view { r -> "RESULT=${r}" }
}
