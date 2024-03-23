process foo {
    input:
    val message

    output:
    val result

    exec:
    result = "$message world"
}

workflow {
    channel.of('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
}