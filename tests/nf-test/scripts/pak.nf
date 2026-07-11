#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process helloPak {
    package "jsonlite", provider: "pak"

    output:
    stdout

    script:
    """
    Rscript -e 'library(jsonlite); cat("Hello pak\\n")'
    """
}

workflow {
    helloPak() | view
}
