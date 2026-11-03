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

    foo.out[0].view { "str: $it" }
    foo.out[1].view { "exp: $it" }
    foo.out[2].view { "x: $it" }
    foo.out[3].view { "y: $it" }
}
