workflow {
  foo()
}

process foo {
  container 'quay.io/nextflow/bash'
  publishDir "gs://rnaseq-nf/scratch/tests", overwrite: true
  output:
  path 'hello.txt'

  script:
  """
  touch hello.txt
  """
}
