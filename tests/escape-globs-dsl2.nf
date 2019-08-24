#!/usr/bin/env nextflow
nextflow.preview.dsl=2
 
process foo {

  input:
  path x

  output:
  path x
  path 'file-\\*.txt'
  path 'file-?.txt' glob false

  '''
  touch file-\\*.txt
  touch file-\\?.txt
  '''

}

workflow {
    Channel.fromPath("$baseDir/data/file\\[a-b\\].txt") | foo

    foo.out[0].println { "match: ${it.name}" }
    foo.out[1].println { "match: ${it.name}" }
    foo.out[2].println { "match: ${it.name}" }
}

