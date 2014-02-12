	Channel
		.from( 4, 1, 7, 5 )
		.sum { it * it } 
		.subscribe {  println "Square: $it" } 
