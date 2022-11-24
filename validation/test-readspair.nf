workflow {
    reads_pair()
    reads_pair
        .out.view()
        .collect()
        .subscribe { assert it.name == ['test.R1.fastq','test.R2.fastq'] }
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

