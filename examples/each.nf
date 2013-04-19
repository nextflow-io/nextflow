#!/usr/bin/env nextflow

channel = channel( 1, 2, 3 )

channel.each {
    println it
}
