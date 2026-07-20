source = channel.of( 'alpha', 'beta', 'delta' )

source.subscribe { str ->
    println "Got: ${str}; len: ${str.length()}"
}