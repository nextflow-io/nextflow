Channel.of( 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 )
    .buffer( size: 3, skip: 2 )
    .view()