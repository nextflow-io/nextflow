process stage1 {
    maxForks 2
    input: val x
    output: val "${x}_s1"
    script:
    """
    sleep 2
    """
}

process stage2 {
    maxForks 2
    input: val x
    output: val "${x}_s2"
    script:
    """
    sleep 2
    """
}

process stage3 {
    maxForks 2
    input: val x
    output: val "${x}_s3"
    script:
    """
    sleep 2
    """
}

process stage4 {
    maxForks 2
    input: val x
    output: val "${x}_s4"
    script:
    """
    sleep 2
    """
}

process stage5 {
    maxForks 2
    input: val x
    output: stdout
    script:
    """
    sleep 2
    echo "done $x"
    """
}

workflow {
    Channel.of(1..25) | stage1 | stage2 | stage3 | stage4 | stage5
}
