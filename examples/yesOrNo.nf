#!/usr/bin/env nextflow

process yesOrNo {
    echo true

    input:
    val ( [1,2] ) as x

    script:
    if( x == 1 ) {
        'echo YES'
    }
    else {
        'echo NO'
    }

}