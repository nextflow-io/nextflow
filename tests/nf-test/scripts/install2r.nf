#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloInstall2r {
    package "jsonlite", provider: "install2r"

    output:
    stdout

    script:
    """
    Rscript -e 'library(jsonlite); cat("Hello install2r\\n")'
    """
}

workflow {
    helloInstall2r() | view
}
