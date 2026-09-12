#!/usr/bin/env nextflow

// A stand-in for `nf-core/fetchngs`: it "downloads" a FASTQ pair for each
// accession. The tools are replaced by dummy files so that the example
// runs anywhere.

nextflow.enable.types = true

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}

params {
    ids: Path       // file of SRA/ENA accessions, one per line
}

process SRATOOLS_FASTERQDUMP {
    tag "${id}"

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
    echo "@${id}/1" > ${id}_1.fastq
    echo "@${id}/2" > ${id}_2.fastq
    """
}

workflow {
    main:
    ch_ids = channel.fromList( params.ids.readLines().findAll { line -> line } )
    ch_samples = SRATOOLS_FASTERQDUMP( ch_ids )

    publish:
    samples = ch_samples
}

output {
    samples: Channel<Sample> { path 'fastq' }
}
