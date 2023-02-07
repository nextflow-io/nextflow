process foo {
  container 'docker.io/pditommaso/my-secret-container:latest'
  debug true

  """
  my-secret-script.sh
  """
}

workflow { 
  foo()
}
