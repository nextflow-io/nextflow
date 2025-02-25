// define a channel emitting three values
source = Channel.of( 'alpha', 'beta', 'delta' )

// subscribe to the channel with a function that prints each value
source.subscribe { str -> println "Got: ${str}; len: ${str.length()}" }