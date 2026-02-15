process upstream {
    maxForks 2
    input: val x
    output: val "${x}_processed"
    script:
    """
    sleep 2
    """
}

process middle {
    maxForks 2
    input: val x
    output: val "${x}_mid"
    script:
    """
    sleep 2
    """
}

process downstream {
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
    Channel.of(1..25) | upstream | middle | downstream
}
