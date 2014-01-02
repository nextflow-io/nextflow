lines = []
lines << 'hello\n' << 'hi\n' << 'hola\n'

result = file('result.txt')

process saveLines {
    echo true

    input:
    file lines as 'item'

    share:
    file result as 'result' to done


    '''
    cat item >> result
    '''

}

done.subscribe { println "Result:\n${it.text}" }
