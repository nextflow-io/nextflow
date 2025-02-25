// mapping function
Channel.of( "hello", "hi", "hey" )
    .min { v -> v.length() }
    .view()