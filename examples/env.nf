#!/usr/bin/env nextflow

echo true
config.env [ 'HELLO' ]  = '1'

task {
    "env | sort"
}

