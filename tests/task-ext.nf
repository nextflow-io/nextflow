
process HELLO {
    debug true

    script:
    """
    echo "args = ${task.ext.args}"
    """
}

workflow {
    HELLO()
}
