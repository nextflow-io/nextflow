#!/usr/bin/env nextflow
// parser_v1: unquoted env inputs/outputs
// parser_v1: implicit process script section

process foo {
    output:
    env FOO
    /FOO=Hello/
}

process bar {
    debug true
    input:
    env FOO
    'echo "bar says $FOO"'
}


workflow {
  foo | bar
}
