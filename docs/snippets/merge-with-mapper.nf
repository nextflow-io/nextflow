odds  = channel.of(1, 3, 5, 7, 9)
evens = channel.of(2, 4, 6)

odds
    .merge( evens ) { a, b -> tuple(b*b, a) }
    .view()