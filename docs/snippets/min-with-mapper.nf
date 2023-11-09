// mapping function
Channel.of( "hello", "hi", "hey" )
    .min { it.size() }
    .view()