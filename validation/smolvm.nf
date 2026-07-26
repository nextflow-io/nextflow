workflow {
  hello()
}

process hello {
  script:
  """
  echo "Hello from smolvm running \$(cat /etc/os-release | grep PRETTY_NAME)"
  """
}
