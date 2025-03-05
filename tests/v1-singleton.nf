#!/usr/bin/env nextflow
// parser_v1: unquoted file output
// parser_v1: implicit process script section

process foo {
  output:
  file x

  '''
  echo -n Hello > x
  '''
}

process bar {
  input:
  file x
  val y

  """
  cat $x
  echo $y
  """

}

workflow {
  foo()
  bar(foo.out, channel.of(1,2,3))
}
