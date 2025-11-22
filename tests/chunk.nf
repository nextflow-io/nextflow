#!/usr/bin/env nextflow

params {
    input: Path
    chunkSize: Integer = 1
}

process foo {
    debug true

    input:
    stdin()

    script:
    "cat -"
}

workflow {
    channel.of(params.input)
        | splitFasta(by: params.chunkSize)
        | foo
}
