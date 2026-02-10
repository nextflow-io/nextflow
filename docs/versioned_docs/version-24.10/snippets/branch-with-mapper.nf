channel.of(1, 2, 3, 40, 50)
    .branch { v ->
        alpha: v < 10
            return v + 2

        beta: v < 50
            return v - 2

        other: true
            return 0
    }
    .set { result }

result.alpha.view { v -> "$v is alpha" }
result.beta.view { v -> "$v is beta" }
result.other.view { v -> "$v is other" }