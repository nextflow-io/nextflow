    Channel
        .from(1,2,3,1,2,3,1)
        .collate( 3, false )
	.subscribe { println it  }

