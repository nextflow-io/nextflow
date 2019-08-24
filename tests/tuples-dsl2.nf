#!/usr/bin/env nextflow
nextflow.preview.dsl=2


process touch {

  input:
    tuple ( id, fileName )
  output:
    tuple ( id, 'file*' )


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    tuple( id, 'file_x' ) 

  output:
    tuple( id, '*') mode flatten

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
        | subscribe { println it }
}
