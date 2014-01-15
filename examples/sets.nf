#!/usr/bin/env nextflow

x = Channel.from( ['a', 'file1'], ['b','file2'] )

process touch {

  input:
    set ( id, fileName ) from x
  output:
    set ( id, 'file*' ) into z


  /
  echo Creating $id
  touch $fileName
  /
}

process makeFiles {
  input:
    set( id, 'file_x' ) from z

  output:
    set( id, '*') into q mode flatten

  /
   cp file_x copy_$id
   touch beta_$id
  /

}

q.subscribe { println it }