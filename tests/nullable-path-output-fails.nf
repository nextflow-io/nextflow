#!/usr/bin/env nextflow

process foo {
  input:
    val id
  output:
    path('output.txt')
  exec:
    println 'hi'
}

workflow {
    foo('foo')
}