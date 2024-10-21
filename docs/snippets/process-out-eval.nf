process sayHello {
    output:
    eval('bash --version')

    script:
    """
    echo Hello world!
    """
}

workflow {
    sayHello | view
}
