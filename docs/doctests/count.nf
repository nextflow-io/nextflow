Channel
    .from(4,1,7,1,1)
    .count(1)
    .subscribe { println "literal $it" }
 // -> 3

Channel
    .from('a','c','c','q','b')
    .count ( ~/c/ )
    .subscribe { println "regexp: $it" }
// -> 2

Channel
    .from('a','c','c','q','b')
    .count { it <= 'c' }
    .subscribe { println "predicate: $it" }
