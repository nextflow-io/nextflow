#!/bin/env nextflow

sequences = queue()

stdin.chunkFasta { str ->
    sequences << it
}

task {
    echo true
    input '-': sequences

    "cat -"
}