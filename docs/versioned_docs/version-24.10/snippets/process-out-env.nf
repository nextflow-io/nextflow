process seq {
    output:
    env 'RESULT'

    script:
    '''
    RESULT=$(seq 5)
    '''
}

workflow {
    seq | view
}
