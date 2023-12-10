process sayHello {
    output:
    cmd('bash --version')

    """
    echo Hello world!
    """
}

workflow {
    sayHello | view
}