Channel.of( 4, 1, 7, 5 )
    .sum { it * it }
    .view { "Square: $it" }