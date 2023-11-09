Channel.of( 1, 2, 3, 4, 5 )
    .filter { it % 2 == 1 }
    .view()