Channel
    .of( 'alpha', 'beta', 'lambda' )
    .subscribe { str ->
        println "Got: ${str}; len: ${str.length()}"
    }