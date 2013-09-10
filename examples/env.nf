#!/usr/bin/env nextflow

/*
 * Shows how manipulate the script execution environment
 */


config.env [ 'HELLO_1' ]  = '1'

process printEnv {
    echo true

    input:
    env 'HELLO_2' using '2'
    env 'HELLO_X' using (['a','b','c'])

    "env | grep HELLO | sort"
}



