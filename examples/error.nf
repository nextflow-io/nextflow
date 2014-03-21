#!/usr/bin/env nextflow

process task1 {
    maxForks 4
    errorStrategy 'ignore'

    input:
    val x from (1,2,3)

    script:
    "echo $x; exit 1"
}

process task2 {
    maxForks 4

    input:
    val x from([4,5,6])

    script:
    "echo $x"

 }
