#!/usr/bin/env nextflow

process foo {
  input:
    val id
  output:
    tuple val(id), path('output.txt')
  exec:
    println 'hi'
}

workflow {
    foo('foo')
}