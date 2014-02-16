    r1 = Channel.create()
    r2 = Channel.create()
    r3 = Channel.create()
	
    Channel
        .from('hello','ciao','hola', 'hi', 'x', 'bonjour')
        .route ( b: r1, c: r2, h: r3 ) { it[0] }

	r3.subscribe { println it }     
