#!/usr/bin/env nextflow

process foo {
  input:
    val id
  output:
    tuple val(id), path('output.txt', nullable: true)
  exec:
    println id
}

process bar {
  input:
    tuple val(id), path(file, nullable: true)
  output:
    val file
  exec:
    sleep 1000L
    println file
}

workflow {
    channel.of('foo') | foo | bar | view()
}