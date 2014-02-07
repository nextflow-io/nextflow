    Channel
        .from( 1, 2, 3, 4, 5 )
        .reduce { a, b ->  println "a: $a	b: $b"; return a+b }
        .subscribe { println "result = $it" }