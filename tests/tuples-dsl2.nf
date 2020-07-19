#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process touch {

  input:
    tuple val(id), val(fileName)
  output:
    tuple val(id), path('file*')


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    tuple val(id), path('file_x')

  output:
    tuple val(id), path('*')

  /
   cp file_x copy_$id
   touch beta_$id
  /
}

workflow {

    Channel
        .from( ['a', 'file1'], ['b','file2'] ) \
        | touch \
        | makeFiles \
        | flatten \
        | subscribe { println it }
}
