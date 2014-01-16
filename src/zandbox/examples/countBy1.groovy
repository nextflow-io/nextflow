import nextflow.Channel

Channel
        .from( 'x', 'y', 'x', 'x', 'z', 'y' )
        .countBy()
        .subscribe { println it }

sleep 100