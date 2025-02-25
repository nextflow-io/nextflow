Channel.of( 1, 2, 3, 1, 2, 3 )
    .buffer { v -> v == 2 }
    .view()