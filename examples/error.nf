#!/usr/bin/env nextflow

process task1 {
    maxForks 4
    errorStrategy 'ignore'

    input:
    val channel(1,2,3) as x

    script:
    "echo $x; exit 1"
}

sleep 500


process task2 {
    maxForks 4

    input:
    val channel(4,5,6) as x

    script:
    "echo $x"

 }