workflow {
  foo()
}

process foo {
  container 'public.cr.seqera.io/mirror/bash'
  publishDir "gs://rnaseq-nf/scratch/tests", overwrite: true
  output:
  path 'hello.txt'

  script:
  """
  touch hello.txt
  """
}
