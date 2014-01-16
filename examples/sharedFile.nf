#!/usr/bin/env nextflow

lines = []
lines << 'hello\n' << 'hi\n' << 'hola\n'

process saveLines {
    echo true

    input:
    file 'item' from lines

    share:
    file result into done


    """
    cat item >> $result
    """

}

done.subscribe {
    println "Result: $it"
    println "\n${it.text}"
}
