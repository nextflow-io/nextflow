#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.chunkSize = 1


process foo {
    echo true

    input:
    stdin()

    "cat -"
}


workflow {
    Channel.from(stdin) \
            | splitFasta( by: params.chunkSize) \
            | foo
}
