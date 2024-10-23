Channel.of( 3, 2, 1, 5, 1, 5 )
    .until { v -> v == 5 }
    .view()