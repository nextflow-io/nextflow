// comparator function
Channel.of( "hello", "hi", "hey" )
    .max { a, b -> a.length() <=> b.length() }
    .view()