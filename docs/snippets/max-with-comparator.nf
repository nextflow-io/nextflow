// comparator function
Channel.of( "hello", "hi", "hey" )
    .max { a, b -> a.size() <=> b.size() }
    .view()