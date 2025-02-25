process cow {
  debug true
  conda 'https://prefix.dev/envs/pditommaso/wave/conda-lock.yml'

  '''
  echo cowpy 'Hello Spack'
  '''
}

workflow {
  cow()
}
