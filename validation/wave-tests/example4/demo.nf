process cow {
  debug true
  conda 'https://prefix.dev/envs/pditommaso/wave/conda-lock.yml'

  script:
  '''
  echo cowpy 'Hello Spack'
  '''
}

workflow {
  cow()
}
