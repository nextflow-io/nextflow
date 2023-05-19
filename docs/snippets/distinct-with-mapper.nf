Channel.of( 1, 1, 2, 2, 2, 3, 1, 1, 2, 4, 6 )
    .distinct { it % 2 }
    .view()