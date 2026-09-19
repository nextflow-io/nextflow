#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloNix {
    package "hello", provider: "nix"

    output:
    stdout

    script:
    """
    hello
    """
}

workflow {
    helloNix() | view
}
