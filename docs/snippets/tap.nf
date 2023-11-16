Channel.of( 'a', 'b', 'c' )
    .tap { log1 }
    .map { it * 2 }
    .tap { log2 }
    .map { it.toUpperCase() }
    .view { "Result: $it" }

log1.view { "Log 1: $it" }
log2.view { "Log 2: $it" }