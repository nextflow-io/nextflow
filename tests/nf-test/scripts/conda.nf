#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloConda {
    package "cowpy", provider: "conda"

    output:
    stdout

    script:
    """
    # cowpy is provided by the conda env; proving it is on PATH confirms activation
    command -v cowpy >/dev/null && echo "Hello conda"
    """
}

workflow {
    helloConda() | view
}
