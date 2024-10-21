// mapping function
Channel.of( "hello", "hi", "hey" )
    .max { v -> v.length() }
    .view()