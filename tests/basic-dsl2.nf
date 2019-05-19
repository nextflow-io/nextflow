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
    file 'input.fa'

    output:
    file 'seq_*'

    """
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    """

}

/* 
 * Simple reverse the sequences 
 */
process reverse {

    input:
    file x 
    
    output:
    stdout()

    """
    cat $x | rev
    """
}

Channel.value(file(params.in)) |
        splitSequences |
        reverse |
        subscribe { println it }

