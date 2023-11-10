Channel
    .of( 'alpha', 'beta', 'lambda' )
    .subscribe { String str ->
        println "Got: ${str}; len: ${str.size()}"
    }