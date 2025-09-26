// Test to validate unstage controls copy error is not printed whe using storeDir. https://github.com/nextflow-io/nextflow/issues/6311

workflow {
  test()
}

process test {
    storeDir "test"

    output:
    path "test.txt"

    script:
    """
    ls > test.txt
    """
}
