#!/usr/bin/env nextflow

params.in = "$HOME/sample.fa"

sequences = file(params.in)
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process splitSequences {

    input:
    file sequences as 'input.fa'

    output:
    file 'seq_*' to records

    """
    $SPLIT input.fa '%^>%' '/^>/' '{*}' -f seq_
    """

}

process reverse {

    input:
    file records as x
    
    output:
    stdout result

    """
    cat $x | rev
    """
}

result.subscribe { println it }
