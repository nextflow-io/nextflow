#!/usr/bin/env nextflow

process foo {
  output:
    path('output.txt', arity: '0..1')
  script:
    true
}

process bar {
  input:
    path(file)
  script:
    true
}

workflow {
    foo | bar
}