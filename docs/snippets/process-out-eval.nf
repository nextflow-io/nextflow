process sayHello {
    output:
    eval('bash --version')

    """
    echo Hello world!
    """
}

workflow {
    sayHello | view
}
