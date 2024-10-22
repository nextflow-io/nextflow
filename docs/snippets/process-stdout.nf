process sayHello {
    output:
    stdout

    script:
    """
    echo Hello world!
    """
}

workflow {
    sayHello | view { "I say... $it" }
}