#!/usr/bin/env nextflow

int count = 0
stdin.chunkFasta { seq ->

    println "seq ($count): " + seq

}


sleep 100