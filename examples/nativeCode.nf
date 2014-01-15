#!/usr/bin/env nextflow

process nativeCode {

    input:
    val x from 'world'

    output:
    val y into stream

    exec:
    y = "Hello $x"
    println "workDir: $workDir"

}


stream.subscribe { println it }

