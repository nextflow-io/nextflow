// mapping function
Channel.of( "hello", "hi", "hey" )
    .max { it.size() }
    .view()