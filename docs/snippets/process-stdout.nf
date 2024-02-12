process sayHello {
    output:
    stdout

    """
    echo Hello world!
    """
}

workflow {
    sayHello | view { "I say... $it" }
}