#!/usr/bin/env nextflow

process splitLetters {

    output:
    file 'chunk_*' into letters mode flatten

    '''
    echo 'Hola' | split -b 1 - chunk_
    '''
}

letters.subscribe onNext: { println it.text }