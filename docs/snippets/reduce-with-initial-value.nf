Channel.of( 1, 2, 3, 4, 5 )
    .reduce( 'result:' ) { acc, v ->
        println acc
        acc + ' ' + v
    }
    .view { result -> "final $result" }