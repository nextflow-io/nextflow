#!/usr/bin/env nextflow

process nativeCode {


    input:
    val x from 'world'

    output:
    val x into stream

    exec:
    y = "Hello $x"
    println "workDir: $workDir"
    workDir.resolve('file.txt').text = y

}


stream.subscribe { println it }

