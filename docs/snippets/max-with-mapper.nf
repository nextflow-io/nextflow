// mapping function
channel.of( "hello", "hi", "hey" )
    .max { v -> v.length() }
    .view()