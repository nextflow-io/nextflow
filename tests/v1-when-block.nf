#!/usr/bin/env nextflow
// parser_v1: mixing script declarations and statements

items = [0,1,2,3,4]
decode = ['zero','one','two','three','fourth']

workflow {
  channel.fromList(items) | foo
  channel.fromList(items) | bar
}

process foo {
    debug true
    tag "${decode[x]}"

    input:
    val x

    when:
    x >= 3

    script:
    """
    echo Foo $x
    """
}

process bar {
    debug true
    tag "${decode[x]}"

    input:
    val x

    when:
    x < 3

    script:
    """
    echo Bar $x
    """
}
