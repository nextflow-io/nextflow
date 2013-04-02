#!/bin/env nextflow

fileName = "${HOME}/Downloads/ampa/multi.fa"
fastaFile = new File(fileName)

seq = queue()
fastaFile.chunkFasta { seq << it }


task ('ampa') {

    //  defines the Input and Output
    input '-':seq
    output result

    // The BASH script to be executed - for each - sequence
    """
    cat - > input.file && AMPA.pl -in=input.file -noplot -rf=result -df=data
    """

}


merge {
    echo true
    input result

    """
    cat $result
    """
}


