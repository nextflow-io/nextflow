
process my_process {
    publishDir "s3://nf-kms-xyz/work/ci-test/publish-s3"

    input:
    val(param)

    output:
    file("HELLO.tsv")

    script:
    """
    echo "Hello, world" > HELLO.tsv
    """
}

workflow {
  Channel.of(1) | my_process
}
