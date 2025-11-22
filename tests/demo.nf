#!/usr/bin/env nextflow

process sayHello {
  debug true
  input:
  val x

  script:
  """
  echo '$x world!'
  """
}

workflow {
  channel.of('Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου') | sayHello
}
