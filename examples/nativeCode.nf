#!/usr/bin/env nextflow

process nativeCode {

    input:
    val 'world' as x

    output:
    val y to stream

    exec:
    y = "Hello $x"

}


stream.subscribe { println it }

