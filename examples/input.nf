#!/usr/bin/env nextflow

x = 1
y = ['a','b']

process foo {
    echo true

    input:
    val x
    val y

    output:
    val y into channel

    "echo $x - $y"

}


channel.subscribe {
    println "foo out: $it"
}
