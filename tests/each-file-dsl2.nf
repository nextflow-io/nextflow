#!/usr/bin/env nextflow
nextflow.preview.dsl=2

process foo {
  echo true

  tag "$x"

  input:
  each file(x)

  """
  grep '>' $x
  """

}


Channel.fromPath("$baseDir/data/p{1,2,3}.fa") | foo