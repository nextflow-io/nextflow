#!/usr/bin/env nextflow

/*
 * Integration test for the unified package management system
 * This test demonstrates the new package directive with backward compatibility
 */

nextflow.enable.dsl=2

// Old style conda directive - should show deprecation warning when preview.package is enabled
process oldStyleConda {
    conda 'samtools=1.17'

    output:
    stdout

    script:
    """
    echo "Old style conda: \$(samtools --version 2>/dev/null | head -1 || echo 'not available')"
    """
}

// Old style pixi directive - should show deprecation warning when preview.package is enabled
process oldStylePixi {
    pixi 'samtools'

    output:
    stdout

    script:
    """
    echo "Old style pixi: \$(samtools --version 2>/dev/null | head -1 || echo 'not available')"
    """
}

// New style package directive with explicit provider
process newStyleExplicit {
    package "samtools=1.17", provider: "conda"

    output:
    stdout

    script:
    """
    echo "New style explicit: \$(samtools --version 2>/dev/null | head -1 || echo 'not available')"
    """
}

// New style package directive with default provider (from config)
process newStyleDefault {
    package "samtools=1.17"

    output:
    stdout

    script:
    """
    echo "New style default: \$(samtools --version 2>/dev/null | head -1 || echo 'not available')"
    """
}

// New style package directive with multiple packages
process newStyleMultiple {
    package ["samtools=1.17", "bcftools=1.18"], provider: "conda"

    output:
    stdout

    script:
    """
    echo "New style multiple: samtools \$(samtools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'n/a'), bcftools \$(bcftools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'n/a')"
    """
}

workflow {
    oldStyleConda() | view
    oldStylePixi() | view
    newStyleExplicit() | view
    newStyleDefault() | view
    newStyleMultiple() | view
}