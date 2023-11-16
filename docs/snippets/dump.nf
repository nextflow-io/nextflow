Channel.of( 1, 2, 3 )
    .map { it+1 }
    .dump(tag: 'foo')

Channel.of( 1, 2, 3 )
    .map { it^2 }
    .dump(tag: 'bar')