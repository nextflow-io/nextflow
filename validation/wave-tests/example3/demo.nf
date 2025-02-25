process cow {
  debug true
  conda 'cowpy=1.1.5'

  '''
  echo cowpy 'Hello Spack'
  '''
}

workflow {
  cow()
}
