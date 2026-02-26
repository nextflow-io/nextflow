// comparator function
channel.of( "hello", "hi", "hey" )
    .max { a, b -> a.length() <=> b.length() }
    .view()