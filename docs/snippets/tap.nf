Channel.of( 'a', 'b', 'c' )
    .tap { log1 }
    .map { v -> v * 2 }
    .tap { log2 }
    .map { v -> v.toUpperCase() }
    .view { result -> "Result: $result" }

log1.view { v -> "Log 1: $v" }
log2.view { v -> "Log 2: $v" }