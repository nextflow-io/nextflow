#!/usr/bin/env nextflow

process sayHello {
  debug true
  input:
  val x

  """
  echo '$x world!'
  """
}

workflow {
  Channel.of('Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου') | sayHello
}
