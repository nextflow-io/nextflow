// comparator function
Channel.of( "hello", "hi", "hey" )
    .min { a, b -> a.size() <=> b.size() }
    .view()