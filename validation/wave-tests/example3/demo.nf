process cow {
  debug true
  conda 'cowpy=1.1.5'

  script:
  '''
  echo cowpy 'Hello Spack'
  '''
}

workflow {
  cow()
}
