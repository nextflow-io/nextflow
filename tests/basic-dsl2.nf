#!/usr/bin/env nextflow
nextflow.preview.dsl=2
 
/* 
 * Command line input parameter 
 */
params.in = "$baseDir/data/sample.fa"


/* 
 * split a fasta file in multiple files 
 */
process splitSequences {

    input:
    path 'input.fa'

    output:
    path 'seq_*'

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """

}

/* 
 * Simple reverse the sequences 
 */
process reverse {

    input:
    path x
    
    output:
    stdout()

    """
    cat $x | rev
    """
}

Channel.value(params.in) |
        splitSequences |
        reverse |
        subscribe { println it }

