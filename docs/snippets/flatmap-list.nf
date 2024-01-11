Channel.of( 1, 2, 3 )
    .flatMap { n -> [ n, n*2, n*3 ] }
    .view()