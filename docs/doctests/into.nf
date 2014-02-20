
     alpha = Channel.create()
     delta = Channel.create()

     Channel
        .from([1,2], ['a','b'], ['p','q'])
        .into( alpha, delta )

     alpha.subscribe { println "first : $it" }
     delta.subscribe { println "second: $it" }
