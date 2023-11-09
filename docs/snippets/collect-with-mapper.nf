Channel.of( 'hello', 'ciao', 'bonjour' )
    .collect { it.length() }
    .view()