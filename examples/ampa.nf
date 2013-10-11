#!/usr/bin/env nextflow

params.in = "${HOME}/sample.fa"

fastaFile = file(params.in)
seq = channel()

/*
 * Splits the input file in chunks containing a single sequences,
 * and send each of it over the 'seq' channel
 */
fastaFile.chunkFasta { seq << it }

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

/*
 * print out each 'result' produced by the above step
 */
result.each {
    println it.text
}
