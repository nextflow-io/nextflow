
process splitLetters {

    output:
    file 'chunk_*' using letters

    '''
    echo 'Hola' | split -b 1 - chunk_
    '''
}

letters.subscribe onNext: { println it.text }