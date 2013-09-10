#!/usr/bin/env nextflow

x = 1
y = ['a','b']

process foo {
    echo true

    input:
    val x
    val y

    output:
    val x using channel

    "echo $x - $y"

}


channel.each {
    println "got: $it"
}
