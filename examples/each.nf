#!/usr/bin/env nextflow

channel1 = new Channel( 10, 20, 30 )
channel2 = new Channel( 'alpha', 'beta', 'gamma' )


channel1.each { println it; sleep 100 }

channel2.eachWithIndex { value, index ->
    println "$index > $value";
    sleep 100
}
