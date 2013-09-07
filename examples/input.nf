#!/usr/bin/env nextflow

x = 1
y = ['a','b']

process foo {
    echo true

    input val: x, from: x
    input y
    output val: x, into: channel

    "echo $x - $y"

}


channel.each {
    println "got: $it"
}
