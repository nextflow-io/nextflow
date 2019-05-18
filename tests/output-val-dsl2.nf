#!/usr/bin/env nextflow
nextflow.preview.dsl=2

x = 100
y = 200

process foo {
  input:
  file fastq

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

foo('dummy')

foo.out[0].println { "str: $it" }
foo.out[1].println { "exp: $it" }
foo.out[2].println { "x: $it" }
foo.out[3].println { "y: $it" }
