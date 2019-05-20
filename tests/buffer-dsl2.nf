#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process blastThemAll {
    echo true

    input:
    file x 

    """
    echo $x
    """
}


Channel
    .fromPath("$baseDir/data/p?.fa") |
    toSortedList |
    flatten |
    buffer(size:2, remainder: true) |
    blastThemAll
