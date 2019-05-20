#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process sayHello {
  echo true
  input:
  val x

  """
  echo '$x world!'
  """
}


Channel.from('Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου') | sayHello
