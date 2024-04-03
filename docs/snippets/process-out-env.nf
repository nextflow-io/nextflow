process myTask {
    output:
    env FOO

    script:
    '''
    FOO=$(ls -a)
    '''
}

workflow {
    myTask | view
}