#!/usr/bin/env nextflow


str = Channel.from('hello', 'hola', 'bonjour', 'ciao')

process printEnv {
    echo true

    input:
    env HELLO from str

    '''
    echo $HELLO world!
    '''
}
