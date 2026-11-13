
process foo {
  input:
  val(index)

  output:
  stdout emit: versions, topic: versions

  script:
  """
  echo 'foo: 0.1.0'
  """
}

process bar {
  input:
  val(index)

  output:
  stdout emit: versions, topic: versions

  script:
  """
  echo 'bar: 0.9.0'
  """
}

workflow {
  channel.fromList( 1..3 ) | foo
  channel.fromList( 1..3 ) | bar

  channel.topic('versions')
  | unique
  | collectFile(name: 'versions.txt', sort: true, storeDir: '.')
}
