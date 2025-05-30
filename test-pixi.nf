#!/usr/bin/env nextflow

process test_pixi {
    pixi 'bwa=0.7.17'
    
    """
    echo "Testing Pixi integration"
    echo "PATH: \$PATH"
    which bwa || echo "bwa not found"
    """
}

workflow {
    test_pixi()
}
