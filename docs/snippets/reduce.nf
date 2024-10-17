Channel.of( 1, 2, 3, 4, 5 )
    .reduce { a, b ->
        println "a: $a b: $b"
        a + b
    }
    .view { result -> "result = $result" }