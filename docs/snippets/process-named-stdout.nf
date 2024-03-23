process sayHello {
    input:
    val cheers

    output:
    stdout emit: verbiage

    script:
    """
    echo -n $cheers
    """
}

workflow {
    things = channel.of('Hello world!', 'Yo, dude!', 'Duck!')
    sayHello(things)
    sayHello.out.verbiage.view()
}