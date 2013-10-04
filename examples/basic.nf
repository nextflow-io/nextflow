#!/usr/bin/env nextflow

params.in = "$HOME/sample.fa"

sequences = file(params.in)
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process splitSequences {

    input:
    file 'input.fa' using sequences

    output:
    file 'seq_*' using records

    """
    $SPLIT input.fa '%^>%' '/^>/' '{*}' -f seq_
    """

}

process reverse {

    input:
        val x using records
    output:
        stdout result

    """
    cat $x | rev
    """
}

result.each { println it }
