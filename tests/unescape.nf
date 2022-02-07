#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

    ```
    awk '/^>/{f="seq_"++d} {print > f}' < input.fa
    ```

}

/*
 * Simple reverse the sequences
 */
process sayHello {

    input:
    path x

    output:
    stdout

    """
    echo ```this is \n an unescaped \n\t string```
    """

}


workflow {
    splitSequences(params.in) | sayHello | view
}

