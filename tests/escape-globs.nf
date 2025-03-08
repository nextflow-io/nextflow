#!/usr/bin/env nextflow

process foo {

  input:
  path x

  output:
  path x
  path 'file-\\*.txt'
  path 'file-?.txt', glob: false

  script:
  '''
  touch file-\\*.txt
  touch file-\\?.txt
  '''

}

workflow {
    Channel.fromPath("$baseDir/data/file\\[a-b\\].txt") | foo

    foo.out[0].view { "match: ${it.name}" }
    foo.out[1].view { "match: ${it.name}" }
    foo.out[2].view { "match: ${it.name}" }
}

