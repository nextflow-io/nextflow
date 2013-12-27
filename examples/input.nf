#!/usr/bin/env nextflow

x = 1
y = ['a','b']

process foo {
    echo true

    input:
    val x
    val y

    output:
    val y using channel

    "echo $x - $y"

}


channel.subscribe {
    println "foo out: $it"
}
