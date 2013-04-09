#!/usr/bin/env nextflow

sequences = channel()

stdin.chunkFasta { str ->
    sequences << str
}

task {
    echo true
    input '-': sequences

    "cat -"
}