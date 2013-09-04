#!/usr/bin/env nextflow

/*
 * Shows how manipulate the script execution environment
 */


config.env [ 'HELLO_1' ]  = '1'

task {
    env HELLO_2: '2'
    input env: 'HELLO_X', from: ['a','b','c']
    echo true

    "env | grep HELLO | sort"
}



