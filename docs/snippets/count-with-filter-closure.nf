Channel.of('a', 'c', 'c', 'q', 'b')
    .count { it <= 'c' }
    .view()