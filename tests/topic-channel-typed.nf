
nextflow.preview.types = true

process foo {
  input:
  index: Integer

  topic:
  record(name: 'foo', version: stdout()) >> 'versions'

  script:
  """
  printf '0.1.0'
  """
}

process bar {
  input:
  index: Integer

  topic:
  record(name: 'bar', version: stdout()) >> 'versions'

  script:
  """
  printf '0.9.0'
  """
}

workflow {
  main:
  ch_inputs = channel.of( 1..3 )
  foo( ch_inputs )
  bar( ch_inputs )

  ch_versions = channel.topic('versions')
    .unique()
    .collect()
    .flatMap { rows ->
      rows.toSorted { r -> r.name }
    }

  publish:
  versions = ch_versions
}

output {
  versions {
    index {
      path 'versions.csv'
    }
  }
}
