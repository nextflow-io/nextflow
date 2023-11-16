process foo {
    input:
    val message

    output:
    val result

    exec:
    result = "$message world"
}

process bar {
    input:
    val message

    output:
    val result

    exec:
    result = message.toUpperCase()
}

workflow {
    channel.of('Hello')
      | map { it.reverse() }
      | (foo & bar)
      | mix
      | view
}