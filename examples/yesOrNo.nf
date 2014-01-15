#!/usr/bin/env nextflow

process yesOrNo {
    echo true

    input:
    val x from 1,2

    script:
    if( x == 1 ) {
        'echo YES'
    }
    else {
        'echo NO'
    }

}