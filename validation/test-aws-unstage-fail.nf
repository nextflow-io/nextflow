process test {
  input: 
	val i
  output: 
        file("test${i}")
        file("test_2_${i}")
  script:
  """
  dd if=/dev/urandom of=test${i} bs=1K count=90
  dd if=/dev/urandom of=test_2_${i} bs=1K count=90
  """
}

workflow {
  channel.of(1) | test
}
