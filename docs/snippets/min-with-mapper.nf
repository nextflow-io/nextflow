// mapping function
channel.of( "hello", "hi", "hey" )
    .min { v -> v.length() }
    .view()