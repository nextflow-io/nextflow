#!/usr/bin/env nextflow

process foo {
  input:
  path fastq

  output:
  val 'Hello'
  val "${fastq.baseName}-${x}.out"
  val x
  val y

  script:
  x = 100
  y = 'two hundred'
  """
  echo bar
  """
}

workflow {
    foo("$baseDir/data/prot.fa")

    foo.out[0].view { str -> "str: $str" }
    foo.out[1].view { exp -> "exp: $exp" }
    foo.out[2].view { x -> "x: $x" }
    foo.out[3].view { y -> "y: $y" }
}
