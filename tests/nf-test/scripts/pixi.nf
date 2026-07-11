#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloPixi {
    package "cowpy", provider: "pixi"

    output:
    stdout

    script:
    """
    # cowpy is provided by the pixi env; proving it is on PATH confirms activation
    command -v cowpy >/dev/null && echo "Hello pixi"
    """
}

workflow {
    helloPixi() | view
}
