#!/usr/bin/env nextflow

nextflow.preview.types = true

process foo {
    input:
    id: String

    output:
    file('output.txt', optional: true)

    script:
    """
    echo ${id}
    """
}

process bar {
    input:
    input: Path?

    stage:
    stageAs 'input.txt', input

    output:
    stdout()

    script:
    '''
    [[ -f input.txt ]] && cat input.txt || echo 'empty input'
    '''
}

workflow {
    channel.of('foo') | foo | bar | view
}
