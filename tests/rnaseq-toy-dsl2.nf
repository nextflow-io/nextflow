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
    file genome
     
    output:
    file 'genome.index*'

    """
    bowtie2-build ${genome} genome.index
    """
}

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {     
    input:
    file genome
    file index
    set pair_id, file(reads)
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam"
 
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
    set pair_id, bam_file
     
    output:
    set pair_id, 'transcripts.gtf'

    """
    cufflinks ${bam_file}
    """
}

/*
 * main flow
 */
genome_file = file(params.genome)
read_pairs = Channel.fromFilePairs( params.reads, checkIfExists: true )

/*
 * main flow
 */
buildIndex(genome_file)
mapping(genome_file, buildIndex.out, read_pairs)
makeTranscript(mapping.out)