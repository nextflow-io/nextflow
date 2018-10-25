#!/usr/bin/env nextflow
/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

/*
 * The reference genome file
 */
genome_file = file(params.genome)

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing 
 * three elements: the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )                                              
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }   
    .set { read_pairs } 
 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    input:
    file genome from genome_file
     
    output:
    file 'genome.index*' into genome_index

    """
    bowtie2-build ${genome} genome.index
    """
}

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {     
    input:
    file genome from genome_file
    file index from genome_index
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam_files
 
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
    set pair_id, bam_file from bam_files
     
    output:
    set pair_id, 'transcripts.gtf' into transcripts

    """
    cufflinks ${bam_file}
    """
}
