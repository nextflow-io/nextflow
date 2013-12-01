import nextflow.Channel

Channel
        .from( 'a', 'b', 'aa', 'bc', 3, 4.5 )
        .grep { it.toString().size() == 1 }
        .subscribe { println it }

sleep 100