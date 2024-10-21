process foo {
  container 'docker.io/pditommaso/my-secret-container:latest'
  debug true

  script:
  """
  my-secret-script.sh
  """
}

workflow { 
  foo()
}
