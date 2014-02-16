    queue1 = Channel.create()
    queue2 = Channel.create()

    Channel
        .from ( 2,4,8 ) 
        .separate( queue1, queue2 ) { a -> [a+1, a*a] }

    queue1.subscribe { println "Channel 1: $it" }
    queue2.subscribe { println "Channel 2: $it" }
