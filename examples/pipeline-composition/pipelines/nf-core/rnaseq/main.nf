#!/usr/bin/env nextflow

// A stand-in for `nf-core/rnaseq`: it aligns each sample, quantifies it, and
// summarizes the run. The tools are replaced by dummy files so that the
// example runs anywhere.

nextflow.enable.types = true

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
    strandedness: String
}

record Alignment {
    id: String
    bam: Path
}

params {
    input: Channel<Sample>              // samplesheet
    aligner: String = 'star_salmon'
    fasta: Path
}

process STAR_ALIGN {
    tag "${sample.id}"

    input:
    sample: Sample
    fasta: Path
    aligner: String

    output:
    record(
        id: sample.id,
        bam: file('*.bam')
    )

    script:
    """
    echo "aligned ${sample.id} against ${fasta.name} with ${aligner} (${sample.strandedness})" > ${sample.id}.bam
    """
}

process SALMON_QUANT {
    tag "${alignment.id}"

    input:
    alignment: Alignment

    output:
    record(
        id: alignment.id,
        counts: file('*.counts.tsv')
    )

    script:
    """
    printf 'gene\\tcount\\nENSG0001\\t42\\n' > ${alignment.id}.counts.tsv
    """
}

process MULTIQC {
    input:
    counts: Bag<Path>

    output:
    file('multiqc_report.html')

    script:
    """
    echo "summarized ${counts.size()} samples" > multiqc_report.html
    """
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
