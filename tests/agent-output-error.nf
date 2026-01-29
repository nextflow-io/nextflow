#!/usr/bin/env nextflow
workflow {
    FAIL()
}


/*
 * Test for agent output mode error handling (NXF_AGENT=true)
 * 
 * This test verifies the error output format in agent mode.
 * Run with: NXF_AGENT=true nextflow run agent-output-error.nf
 */

process FAIL {
    script:
    """
    echo "Error message here" >&2
    exit 127
    """
}
