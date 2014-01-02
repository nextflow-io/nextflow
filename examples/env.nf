#!/usr/bin/env nextflow

/*
 * Shows how manipulate the script execution environment
 */


config.env [ 'HELLO_1' ]  = '1'

process printEnv {
    echo true

    input:
    env '2' as 'HELLO_2'
    env (['a','b','c']) as 'HELLO_X'

    "env | grep HELLO | sort"
}



