#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process blastThemAll {
    echo true

    input:
    path x

    """
    echo $x
    """
}


workflow {
    Channel
        .fromPath("$baseDir/data/p?.fa") |
        toSortedList |
        flatten |
        buffer(size:2, remainder: true) |
        blastThemAll
}
