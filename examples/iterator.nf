#!/usr/bin/env nextflow

list1 = [1,2]
list2 = ['Hola', 'mundo']
list3 = ['alpha','beta','delta']

process hola {
    echo true

    input:
    val list1 as x
    each list2 as y
    each list3 as z

    """
    echo 'x: $x; y: $y; z: $z'
    """

}
