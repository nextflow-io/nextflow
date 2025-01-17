Channel.of( 1, 2, 3, 4, 5 )
    .filter { v -> v % 2 == 1 }
    .view()