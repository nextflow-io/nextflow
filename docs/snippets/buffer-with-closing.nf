Channel.of( 1, 2, 3, 1, 2, 3 )
    .buffer { it == 2 }
    .view()