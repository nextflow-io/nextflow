#!/usr/bin/env nextflow

/*
 * Shows how manipulate the script execution environment
 */


echo true
config.env [ 'HELLO_1' ]  = '1'

task {
    environment HELLO_2: 2

    "env | sort"
}



