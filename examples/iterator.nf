#!/usr/bin/env nextflow

list1 = [1,2]
list2 = ['Hola', 'mundo']
list3 = ['alpha','beta','delta']

process hola {
    echo true

    input:
    val x from list1
    each y from list2
    each z from list3

    """
    echo 'x: $x; y: $y; z: $z'
    """

}
