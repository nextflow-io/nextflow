channel.of( 1, 2, 3 )
    .multiMap { v -> alpha: beta: v }
    .set { result }

result.alpha.view { v -> "alpha $v" }
result.beta.view { v -> "beta $v" }