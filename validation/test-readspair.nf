workflow {
  reads_pair()
  reads_pair.out
    .view()
    .subscribe { files ->
      assert files*.name == ['test.R1.fastq','test.R2.fastq']
    }
}

process reads_pair {
  output:
    file("reads/*")

  script:
    """
    mkdir reads
    touch reads/test.R1.fastq reads/test.R2.fastq
    """
}

