
process foo {
    input:
    val n_files

    output:
    path '*.fastq'

    script:
    """
    for i in `seq 1 ${n_files}`; do
        touch sample_\${i}.fastq
    done
    """
}

process bar {
    input:
    path '*.fastq'

    output:
    stdout

    script:
    """
    for f in `ls *.fastq`; do
        echo \$f
        cat \$f
    done
    """
}

workflow {
    foo(10) | bar
}
