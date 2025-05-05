// comparator function
channel.of( "hello", "hi", "hey" )
    .min { a, b -> a.length() <=> b.length() }
    .view()