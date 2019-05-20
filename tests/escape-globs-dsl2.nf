#!/usr/bin/env nextflow
nextflow.preview.dsl=2
 
process foo {

  input:
  file x

  output:
  file x
  file 'file-\\*.txt'
  file 'file-?.txt' glob false 

  '''
  touch file-\\*.txt
  touch file-\\?.txt
  '''

}

Channel.fromPath("$baseDir/data/file\\[a-b\\].txt") | foo

foo.out[0].println { "match: ${it.name}" }
foo.out[1].println { "match: ${it.name}" }
foo.out[2].println { "match: ${it.name}" }
