#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloUv {
    package "cowsay", provider: "uv"

    output:
    stdout

    script:
    """
    python -c "import cowsay; print('Hello uv')"
    """
}

workflow {
    helloUv() | view
}
