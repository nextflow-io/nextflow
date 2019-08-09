#!/usr/bin/env nextflow
nextflow.preview.dsl=2

x = 100
y = 200

process foo {
  input:
  path fastq

  output:
  val 'Hello'
  val "${fastq.baseName}-${x}.out"
  val x
  val y

  script:
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
