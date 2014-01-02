#!/usr/bin/env nextflow

params.in = '~/sample.fa'
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process split {
    input:
    file  file(params.in) as 'query.fa'

    output:
    file 'seq_*' to splits

    """
    $SPLIT query.fa '%^>%' '/^>/' '{*}' -f seq_
    """
}


process printTwo {
    echo true

    input:
    file splits as 'chunk'

    output:
    file 'chunk1:chunk3' to two_chunks flat true

    """
    cat chunk* | rev
    """

}

process printLast {
    echo true

    input:
    file two_chunks as 'chunk'

    output:
    file 'chunk' to result

    """
    cat chunk
    """
}
