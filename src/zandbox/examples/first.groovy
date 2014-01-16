import nextflow.Channel

Channel.from( 1, 2, 3 ).first() .subscribe { println it }

//
Channel.from( 1, 2, 'a', 'b', 3 ) .first( String ) .subscribe { println it }

//
Channel.from( 'a', 'aaa', 'aaa' ) .first( ~/aa.*/) .subscribe { println it }

Channel.from( 1,2,3,4,5 ) .first { it > 3 } .subscribe { println it }


sleep 100