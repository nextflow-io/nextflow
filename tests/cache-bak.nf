workflow {
  foo()
}

process foo {
  debug true
  script:
  /echo Hello world/
}