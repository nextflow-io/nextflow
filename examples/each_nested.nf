#!/usr/bin/env nextflow

channel = channel( 1, 2, 3 )
channel2 = channel( 4, 5, 6 )


channel.each { it1 ->

    channel2.each { it2 -> println "$it1 - $it2" }

}


sleep 1000