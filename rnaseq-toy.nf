#!/usr/bin/env nextflow
/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
