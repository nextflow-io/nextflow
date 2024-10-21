#!/usr/bin/env nextflow

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

    script:
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
    stdout

    script:
    """
    cat $x | rev
    """
}


workflow {
    splitSequences(params.in) | reverse | view
}

