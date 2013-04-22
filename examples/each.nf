#!/usr/bin/env nextflow

channel = new Channel( 1, 2, 3 )

channel.each {
    println it
}
