#!/usr/bin/env nextflow

nextflow.preview.types = true

process TOUCH {
    input:
    id: String

    output:
    record(
        id: id,
        fastq_1: file('*_1.fastq'),
        fastq_2: file('*_2.fastq')
    )

    script:
    """
    touch ${id}_1.fastq
    touch ${id}_2.fastq
    """
}

process FASTQC {
    input:
    sample: Sample

    output:
    record(
        id: sample.id,
        html: file('*.html'),
        zip: file('*.zip')
    )

    script:
    """
    touch ${sample.id}.html
    touch ${sample.id}.zip
    """
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}

workflow {

    ch_samples = TOUCH( channel.of('a', 'b', 'c') )
    ch_fastqc = FASTQC(ch_samples)
    ch_fastqc.view()
}
