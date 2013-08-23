#!/usr/bin/env nextflow

x = 1
y = ['a','b']

task {
    echo true
    input x
    input y

    "echo $x - $y"

}
