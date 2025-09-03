#!/usr/bin/env nextflow

/*
 * Test script for the unified package management system
 */

nextflow.enable.dsl=2

// Test the new package directive with conda provider
process testConda {
    package "samtools=1.17", provider: "conda"

    output:
    stdout

    script:
    """
    samtools --version | head -1
    """
}

// Test the new package directive with pixi provider
process testPixi {
    package "samtools=1.17", provider: "pixi"

    output:
    stdout

    script:
    """
    samtools --version | head -1
    """
}

// Test the new package directive with default provider
process testDefault {
    package "samtools=1.17"

    output:
    stdout

    script:
    """
    samtools --version | head -1
    """
}

workflow {
    testConda() | view { "Conda: ${it.trim()}" }
    testPixi() | view { "Pixi: ${it.trim()}" }
    testDefault() | view { "Default: ${it.trim()}" }
}