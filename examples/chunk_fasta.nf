#!/bin/env nextflow

sequences = stdin.chunkFasta()

task {
    echo true
    input '-': sequences

    "cat -"
}