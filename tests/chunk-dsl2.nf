#!/usr/bin/env nextflow

params.chunkSize = 1

process foo {
    debug true

    input:
    stdin()

    script:
    "cat -"
}

workflow {
    Channel.from(stdin) \
            | splitFasta( by: params.chunkSize) \
            | foo
}
