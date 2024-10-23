// comparator function
Channel.of( "hello", "hi", "hey" )
    .min { a, b -> a.length() <=> b.length() }
    .view()