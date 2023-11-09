Channel.of( 1, 2, 3, 4, 5 )
    .reduce( 'result:' ) { accum, v ->
        println accum
        accum + ' ' + v
    }
    .view { "final $it" }