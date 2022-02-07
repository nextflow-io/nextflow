#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.in = "$baseDir/data/sample.fa"

include { sayHello } from './unescape'

workflow {
    sayHello(params.in) | view
}