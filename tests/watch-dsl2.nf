#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.events = 'create'
params.files = 'examples/data/*.fa'


process align {
  input:
  file fasta

  output:
  file aln

  """
  t_coffee -in $fasta 1> aln
  """
}

/*
 * main flow
 */
 
Channel
    .watchPath(params.files, params.events) \
    | align \
    | subscribe {
          println '------'
          println it.text
        }
