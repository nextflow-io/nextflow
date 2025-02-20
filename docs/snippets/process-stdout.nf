process sayHello {
    output:
    stdout

    script:
    """
    echo Hello world!
    """
}

workflow {
    sayHello | view { message -> "I say... $message" }
}