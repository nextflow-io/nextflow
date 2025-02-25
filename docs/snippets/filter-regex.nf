Channel.of( 'a', 'b', 'aa', 'bc', 3, 4.5 )
    .filter( ~/^a.*/ )
    .view()