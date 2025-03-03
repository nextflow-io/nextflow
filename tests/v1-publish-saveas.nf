#!/usr/bin/env nextflow
// parser_v1: implicit environment variables
// parser_v1: method pointers
// parser_v1: implicit process script section

def rule( file ) {
  if( file == 'file_1.txt' )
    return "alpha/$file"

  if( file == 'file_2.txt' )
    return null

  if( file == 'file_3.txt' )
     return "$PWD/results/gamma/$file"

}

process foo {
  publishDir path: 'results', saveAs: this.&rule

  input: each x
  output: path '*.txt'
  """
  touch file_${x}.txt
  """

}

workflow {
  foo([1,2,3])
}
