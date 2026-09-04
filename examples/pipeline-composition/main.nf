#!/usr/bin/env nextflow

// A meta-pipeline: it fetches FASTQ samples with `nf-core/fetchngs` and
// analyzes them with `nf-core/rnaseq`, composing the two pipelines with
// regular dataflow logic.
//
// The params and output blocks of each included pipeline are imported as
// record types, so that the meta-pipeline declares one param per pipeline
// instead of replicating every param.

nextflow.enable.types = true

include {
    params as FetchngsParams ;
    workflow as NFCORE_FETCHNGS
} from './pipelines/nf-core/fetchngs'

include {
    params as RnaseqParams ;
    workflow as NFCORE_RNASEQ ;
    output as RnaseqOutput
} from './pipelines/nf-core/rnaseq'

params {
    fetchngs: FetchngsParams        // ids
    strandedness: String = 'auto'   // unique to the meta-pipeline
    rnaseq: RnaseqParams            // input, aligner, fasta
}

workflow {
    main:
    // fetch FASTQ samples from NCBI SRA
    fetchngs = NFCORE_FETCHNGS( params.fetchngs )

    // adapt fetchngs output to rnaseq input (add strandedness)
    ch_samples = fetchngs.samples.map { sample ->
        sample + record(strandedness: params.strandedness)
    }

    // perform RNAseq analysis (ch_samples overrides params.rnaseq.input)
    rnaseq = NFCORE_RNASEQ( params.rnaseq + record(input: ch_samples) )

    publish:
    rnaseq = rnaseq
}

output {
    rnaseq: RnaseqOutput {}
}
