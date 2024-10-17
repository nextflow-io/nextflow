Channel.of( 4, 1, 7, 5 )
    .sum { v -> v * v }
    .view { result -> "Square: $result" }