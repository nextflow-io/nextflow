#!/usr/bin/env nextflow

echo true

process sayhello {

    input:
    val (['hello','hi']) as x

    output:
    val y to all

    share:
    val ('world') as y to z

    """
    echo '$x $y!'
    """
}

z.subscribe { println "Complete: $it" }

all.subscribe { println "All: $it" }


process inc {
    echo true

    input:
        val ([1,2,3,4]) as time

    share:
        val 0 as x
        val 10 as y
        val '.' as w


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