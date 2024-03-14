#!/usr/bin/env nextflow

process foo {
  input:
    val bar
    val baz
  output: 
    stdout

  script:
    """
    echo $bar
    echo $baz
    """
}


workflow foo_wrapper {
  take:
    bar
    baz
  main:
    foo(bar: bar, baz: baz)
  emit:
    foo.out
}


workflow {
  foo_wrapper(bar: 'bar', baz: 'baz')
    | view
}