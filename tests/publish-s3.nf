process my_process {
    publishDir "s3://nextflow-ci/work/ci-test/publish-s3"

    input:
    val(param) from Channel.from(1)

    output:
    file("HELLO.tsv") into output_ch

    script:
    """
    echo "Hello, world" > HELLO.tsv
    """
}