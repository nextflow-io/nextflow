channel.of( 1, 2, 3 )
    .map { v -> v + 1 }
    .dump(tag: 'plus1')

channel.of( 1, 2, 3 )
    .map { v -> v ** 2 }
    .dump(tag: 'exp2')