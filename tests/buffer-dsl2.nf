#!/usr/bin/env nextflow

process blastThemAll {
    debug true

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
