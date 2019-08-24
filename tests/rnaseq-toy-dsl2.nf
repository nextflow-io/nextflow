#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"


/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    input:
    path genome
     
    output:
    path 'genome.index*'

    """
    bowtie2-build ${genome} genome.index
    """
}

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {     
    input:
    path genome
    path index
    tuple pair_id, path(reads)
 
    output:
    tuple pair_id, path("tophat_out/accepted_hits.bam")
 
    """
    tophat2 genome.index ${reads}
    """
}

/*
 * Step 3. Assembles the transcript by using the "cufflinks" 
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {
    publishDir "results"
    
    input:
    tuple pair_id, path(bam_file)
     
    output:
    tuple pair_id, path('transcripts.gtf')

    """
    cufflinks ${bam_file}
    """
}

/*
 * main flow
 */
read_pairs = Channel.fromFilePairs( params.reads, checkIfExists: true )

/*
 * main flow
 */
workflow {
    buildIndex(params.genome)
    mapping(params.genome, buildIndex.out, read_pairs)
    makeTranscript(mapping.out)
}
