process reads_pair {

output:
  file("reads/*") into ch_out

script:
  """
  mkdir reads
  touch reads/test.R1.fastq reads/test.R2.fastq
  """
}

assert ch_out.view().val.collect { it.name } == ['test.R1.fastq','test.R2.fastq'] 