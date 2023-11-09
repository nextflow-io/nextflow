Channel.of( 3, 2, 1, 5, 1, 5 )
    .until { it == 5 }
    .view()