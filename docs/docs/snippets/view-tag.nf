channel.of( 1, 2, 3 )
    .map { v -> v + 1 }
    .view(tag: 'plus1')

channel.of( 1, 2, 3 )
    .map { v -> v ** 2 }
    .view(tag: 'exp2')
