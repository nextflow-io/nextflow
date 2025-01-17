Channel.of(1, 2, 3)
    .map { v -> [v, v*v] }
    .view { num, sqr -> "The square of $num is $sqr" }