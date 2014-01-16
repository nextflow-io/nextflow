import nextflow.Channel

Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .grep( Number )
        .subscribe { println it }

sleep 100