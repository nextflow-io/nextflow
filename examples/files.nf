#!/usr/bin/env nextflow

params.in = '~/sample.fa'
SPLIT = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process split {
    input:
    file 'query.fa' from file(params.in)

    output:
    file 'seq_*' into splits

    """
    $SPLIT query.fa '%^>%' '/^>/' '{*}' -f seq_
    """
}


process printTwo {
    echo true

    input:
    file 'chunk' from splits

    output:
    file 'chunk1:chunk3' into two_chunks mode flatten

    """
    cat chunk* | rev
    """

}

process printLast {
    echo true

    input:
    file 'chunk' from two_chunks

    output:
    file 'chunk' into result

    """
    cat chunk
    """
}
