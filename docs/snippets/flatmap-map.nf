Channel.of( 1, 2, 3 )
    .flatMap { n -> [ number: n, square: n*n, cube: n*n*n ] }
    .view { entry -> "${entry.key}: ${entry.value}" }