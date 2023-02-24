#!/usr/bin/env nextflow

process foo {
  input:
    val metadata
    path ('star/*'), defaultValue: []
    path ('hisat2/*'), defaultValue: []
    path ('salmon/*'), defaultValue: []
  output:
    stdout
  script:
    """
    echo 'metadata: ${metadata}'
    [[ -d star ]] && ls star || echo 'skipping star directory'
    [[ -d hisat2 ]] && ls hisat2 || echo 'skipping hisat2 directory'
    [[ -d salmon ]] && ls salmon || echo 'skipping salmon directory'
    """
}

workflow {
    metadata = Channel.of('foo')
    foo(metadata) | view
}