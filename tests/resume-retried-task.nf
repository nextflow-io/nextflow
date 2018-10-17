#!/usr/bin/env nextflow

process foo {
    errorStrategy 'retry'
    maxRetries 3 

    script:
    if( task.attempt < 3 )
	"""
    exit 1 
	"""
    else 
    """
    echo ciao
    """
}
