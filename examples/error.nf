#!/usr/bin/env nextflow

process task1 {
    maxForks 4
    errorStrategy 'ignore'

    input:
    val x using channel(1,2,3)

    script:
    "echo $x; exit 1"
}

sleep 500


process task2 {
    maxForks 4

    input:
    val x using channel(4,5,6)

    script:
    "echo $x"

 }