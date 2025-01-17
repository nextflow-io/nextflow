process sayHello {
    output:
    eval('echo Hello world!')

    script:
    """
    true
    """
}

workflow {
    sayHello | view
}
