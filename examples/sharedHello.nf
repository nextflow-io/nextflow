#!/usr/bin/env nextflow

echo true

process sayhello {

    input:
    val x from (['hello','hi'])

    output:
    val y into all

    share:
    val y from 'world' into z

    """
    echo '$x $y!'
    """
}

z.subscribe { println "Complete: $it" }

all.subscribe { println "All: $it" }


process inc {
    echo true

    input:
        val time from 1,2,3,4

    share:
        val x from 0
        val y from 10
        val w from '.'


    script:
        x++
        y++
        w+='.'

        """
        echo x: $x
        echo y: $y
        echo w: $w
        """

}