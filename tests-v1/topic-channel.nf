
nextflow.preview.topic = true

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
  Channel.of( 1..3 ) | foo
  Channel.of( 1..3 ) | bar

  Channel.topic('versions')
  | unique
  | collectFile(name: 'versions.txt', sort: true, storeDir: '.')
}
