#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Demonstrates an explicit relative path to a manifest file in the module
// directory (the previously-supported behavior), resolved via the uv provider.
process helloUvManifest {
    package "${projectDir}/../fixtures/requirements.txt", provider: "uv"

    output:
    stdout

    script:
    """
    python -c "import cowsay; print('Hello uv (requirements.txt)')"
    """
}

workflow {
    helloUvManifest() | view
}
