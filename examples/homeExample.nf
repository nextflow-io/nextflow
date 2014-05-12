#!/usr/bin/env nextflow

str = Channel.from('hello', 'hola', 'bonjour', 'ciao').map { it+'\n' }

process printAll {
   echo true

   input:
   stdin str

   """
   cat -
   """
}	