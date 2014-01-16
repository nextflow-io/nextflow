#!/usr/bin/env nextflow


str = Channel.from('hello', 'hola', 'bonjour', 'ciao')

process printEnv {
    echo true

    input:
    env str as HELLO

    '''
    echo $HELLO world!
    '''
}
