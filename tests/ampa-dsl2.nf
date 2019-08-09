#!/usr/bin/env nextflow 
nextflow.preview.dsl=2

params.in = "$baseDir/data/sample.fa"


/*
 * For each sequence that is sent over the 'seq' channel
 * the below task is executed
 */
process ampaTask {

    input:
    path seq

    output:
    path 'result'

    // The BASH script to be executed - for each - sequence
    """
    AMPA.pl -in=${seq} -noplot -rf=result -df=data
    """

}

workflow {
    Channel.fromPath(params.in) |
            splitFasta(file:true) |
            ampaTask |
            view { it.text }
}

