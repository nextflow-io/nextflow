process slow {
    maxForks 1
    input: val x
    output: stdout
    script:
    """
    sleep 2
    echo "done $x"
    """
}

workflow {
    Channel.of(1..50) | slow
}
