#!/usr/bin/env nextflow

process foo {
  input:
    val id
  output:
    path('output.txt', nullable: true)
  exec:
    println id
}

process bar {
  input:
    path(file, nullable: true)
  output:
    val file
  exec:
    sleep 1000L
    println file
}

workflow {
    channel.of('foo') | foo | bar | view()
}