channel.of('a', 'c', 'c', 'q', 'b')
    .count { v -> v <= 'c' }
    .view()