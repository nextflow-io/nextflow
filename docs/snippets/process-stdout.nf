process hello {
    output:
    stdout

    script:
    """
    echo "Hello world!"
    """
}

workflow {
    hello | view { message -> "I say... $message" }
}