
    Channel.from ( 1, 2, 3 )
           .flatMap { it -> [ number: it, square: it*it ] }
           .subscribe { println it.key + ': ' + it.value }