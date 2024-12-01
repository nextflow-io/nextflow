process myTask {
    output:
    env 'FOO'

    script:
    '''
    FOO=$(seq 5)
    '''
}

workflow {
    myTask | view
}
