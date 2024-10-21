#!/usr/bin/env nextflow 

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

    script:
    """
    AMPA.pl -in=${seq} -noplot -rf=result -df=data
    """

}

workflow {
    Channel.fromPath(params.in)
        | splitFasta(file:true)
        | ampaTask
        | view { file -> file.text }
}

