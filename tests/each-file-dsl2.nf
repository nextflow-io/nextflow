#!/usr/bin/env nextflow

process foo {
  echo true
  tag "$x"

  input:
  each path(x)

  """
  grep '>' $x
  """
}


workflow {
    Channel.fromPath("$baseDir/data/p{1,2,3}.fa") | foo
}
