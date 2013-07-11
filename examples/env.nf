#!/usr/bin/env nextflow

/*
 * Shows how manipulate the script execution environment
 */


config.env [ 'HELLO_1' ]  = '1'

task {
    env HELLO_2: 2
    echo true

    "env | sort"
}



