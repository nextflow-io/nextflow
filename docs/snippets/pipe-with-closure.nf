process foo {
    input:
    val message
    val suffix

    output:
    val result, emit: suffixed

    exec:
    result = "${message}${suffix}"
}

workflow {
    suffix = ' world!'
    channel.from('Hello','Hola','Ciao')
      | map { it.toUpperCase() }
      | { _ -> foo(_, suffix) }
      | view
}