Channel.of( 'hello', 'ciao', 'bonjour' )
    .collect { v -> v.length() }
    .view()