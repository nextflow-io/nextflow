#!/usr/bin/env nextflow


process foo {
  input:
    val id
  output:
    map( [id: val(id), first: path('first.txt'), second: path('second.txt')], emit: 'records' )
  script:
    """
    echo 'first' > first.txt
    echo 'second' > second.txt 
    """
}

process bar {
  input:
    map( [id: val(id), first: path('first.txt'), second: path('second.txt')] )
  output:
    stdout
  script:
    """
    echo ${id}
    cat first.txt
    cat second.txt
    """
}

workflow {
  Channel.of(1, 2, 3)
    | foo
    | bar
    | view
}
