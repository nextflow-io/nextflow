Channel.of(1, 2, 3)
    .map { it -> [it, it*it] }
    .view { num, sqr -> "The square of $num is $sqr" }