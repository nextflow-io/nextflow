process hello {
    output:
    eval('echo Hello world!')

    script:
    """
    true
    """
}

workflow {
    hello | view
}
