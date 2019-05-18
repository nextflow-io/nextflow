#!/usr/bin/env nextflow 
nextflow.preview.dsl=2

params.in = "$baseDir/data/sample.fa"


/*
 * For each sequence that is sent over the 'seq' channel
 * the below task is executed
 */
process ampaTask {

    input:
    file seq

    output:
    file result

    // The BASH script to be executed - for each - sequence
    """
    AMPA.pl -in=${seq} -noplot -rf=result -df=data
    """

}

Channel.fromPath(params.in) |
            splitFasta |
            ampaTask |
            view { it.text }

