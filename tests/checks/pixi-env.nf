#!/usr/bin/env nextflow
workflow {
    sayHello() | view
}


/*
 * Test for Pixi environment support
 */

process sayHello {
    pixi 'cowpy'

    output:
    stdout

    script:
    """
    cowpy "hello pixi"
    """
}
