#!/usr/bin/env nextflow

/*
 * Show how get data from the 'stdin' virtual file
 *
 * The data is supposed to be a FASTA file which is splitted
 * is chunks containing a single sequence
 *
 * The number of sequences in each each is controlled by the command
 * line parameter '--chunkSize' (--chunk-size is a synonym for the same)
 */


params.chunkSize = 1
sequences = new Channel()

stdin.chunkFasta( params.chunkSize ) { str ->
    sequences << str
}

task {
    echo true
    input '-': sequences

    "cat -"
}