#!/usr/bin/env nextflow

nextflow.enable.types = true

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
    stageAs input, 'input.txt'

    output:
    stdout()

    script:
    '''
    [[ -f input.txt ]] && cat input.txt || echo 'empty input'
    '''
}

workflow {
    bar(foo(channel.of('foo'))).view()
}
