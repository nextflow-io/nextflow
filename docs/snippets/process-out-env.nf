process myTask {
    output:
    env FOO

    script:
    '''
    FOO=$(ls -a1)
    '''
}

workflow {
    myTask | view
}
