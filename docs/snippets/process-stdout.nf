process sayHello {
    output:
    stdout

    """
    echo Hello world!
    """
}

workflow {
    sayHello | view { message -> "I say... $message" }
}