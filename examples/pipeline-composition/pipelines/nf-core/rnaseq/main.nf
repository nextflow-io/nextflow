#!/usr/bin/env nextflow

// A stand-in for `nf-core/rnaseq`: it aligns each sample, quantifies it, and
// summarizes the run. The tools are replaced by dummy files so that the
// example runs anywhere.

nextflow.enable.types = true

include { STAR_ALIGN } from './modules/nf-core/star/align'
include { SALMON_QUANT } from './modules/nf-core/salmon/quant'
include { MULTIQC } from './modules/nf-core/multiqc'

params {
    input: Channel<Sample>              // samplesheet
    aligner: String = 'star_salmon'
    fasta: Path
}

workflow {
    main:
    ch_bams = STAR_ALIGN( params.input, params.fasta, params.aligner )
    ch_counts = SALMON_QUANT( ch_bams )
    val_multiqc = MULTIQC( ch_counts.map { c -> c.counts }.collect() )

    publish:
    bams    = ch_bams.map { a -> a.bam }
    counts  = ch_counts.map { c -> c.counts }
    multiqc = val_multiqc
}

output {
    bams: Channel<Path> { path 'bams' }
    counts: Channel<Path> { path 'counts' }
    multiqc: Path { path 'multiqc' }
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
    strandedness: String
}
