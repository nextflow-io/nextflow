#!/usr/bin/env nextflow
workflow {
    sayHello() | view
}


/*
 * Test for Pixi environment support
 */

process sayHello {
    pixi 'python=3.8'

    output:
    stdout

    script:
    '''
    python -c "import sys; print(f'Hello from Python {sys.version_info.major}.{sys.version_info.minor}!')"
    '''
}
