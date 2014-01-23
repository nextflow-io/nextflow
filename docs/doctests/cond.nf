    x = Channel.from( 'a', 'b', 'c')   

    process simpleSum {
        input: 
        val x 

        exec:
        println "Hello Mr. $x"
    }
