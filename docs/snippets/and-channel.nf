process foo {
    input:
    val begin
    val middle
    val end

    output:
    val result

    exec:
    result = "$begin $middle $end"
}

workflow {
    ch_begin = channel.of('Hello')
    ch_middle = channel.of('world')
    ch_end = channel.of('!!!')

    (ch_begin & ch_middle & ch_end)
      | foo
      | view
}