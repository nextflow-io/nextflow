#!/usr/bin/env nextflow

// A stand-in for `nf-core/fetchngs`: it "downloads" a FASTQ pair for each
// accession. The tools are replaced by dummy files so that the example
// runs anywhere.

nextflow.enable.types = true

include { SRATOOLS_FASTERQDUMP } from './modules/nf-core/sratools/fasterqdump'

params {
    ids: Path       // file of SRA/ENA accessions, one per line
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

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}
