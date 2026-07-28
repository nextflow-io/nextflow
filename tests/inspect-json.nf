#!/usr/bin/env nextflow

process foo {
  container 'ubuntu:latest'

  input:
  val x

  output:
  stdout

  script:
  """
  echo $x
  """
}

workflow {
  channel.of('hello') | foo
}
