#!/usr/bin/env nextflow

lines = []
lines << 'hello\n' << 'hi\n' << 'hola\n'

process saveLines(echo: true) {

    input:
    file lines as 'item'

    share:
    file result  to done


    """
    cat item >> $result
    """

}

done.subscribe {
    println "Result: $it"
    println "\n${it.text}"
}
