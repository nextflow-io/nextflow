odds  = Channel.of(1, 3, 5, 7, 9)
evens = Channel.of(2, 4, 6)

odds
    .merge( evens ) { a, b -> tuple(b*b, a) }
    .view()