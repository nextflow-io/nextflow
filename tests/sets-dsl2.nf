#!/usr/bin/env nextflow
nextflow.preview.dsl=2


process touch {

  input:
    set ( id, fileName )
  output:
    set ( id, 'file*' )


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    set( id, 'file_x' ) 

  output:
    set( id, '*') mode flatten

  /
   cp file_x copy_$id
   touch beta_$id
  /
}


Channel
    .from( ['a', 'file1'], ['b','file2'] ) \
    | touch \
    | makeFiles \
    | subscribe { println it }
