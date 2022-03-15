#!/usr/bin/env nextflow

process sayHello {
  echo true
  input:
  val x

  """
  echo '$x world!'
  """
}


workflow {
    Channel.from('Bojour', 'Ciao', 'Hello', 'Hola', 'Γεια σου') | sayHello
}
