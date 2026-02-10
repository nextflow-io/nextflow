channel.of( 1, 2, 3, 4 )
    .multiMap { v ->
        inc: v + 1
        squared: v * v
    }
    .set { result }

result.inc.view { v -> "inc $v" }
result.squared.view { v -> "squared $v" }