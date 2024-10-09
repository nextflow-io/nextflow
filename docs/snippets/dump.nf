Channel.of( 1, 2, 3 )
    .map { v -> v + 1 }
    .dump(tag: 'foo')

Channel.of( 1, 2, 3 )
    .map { v -> v ^ 2 }
    .dump(tag: 'bar')