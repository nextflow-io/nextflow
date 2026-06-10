#!/usr/bin/env nextflow
workflow {
    HELLO | WORLD | view
}


/*
 * Test for agent output mode (NXF_AGENT_MODE=true)
 * 
 * This test verifies the minimal agent-friendly output format.
 * Run with: NXF_AGENT_MODE=true nextflow run agent-output.nf
 */

process HELLO {
    output:
    stdout

    script:
    """
    echo "Hello from agent mode"
    """
}

process WORLD {
    input:
    val msg

    output:
    stdout

    script:
    """
    echo "World received: ${msg}"
    """
}
