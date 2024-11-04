Channel.of(1, 2, 3, 40, 50)
    .branch { v ->
        small: v < 10
        large: v < 50
        other: true
    }
    .set { result }

result.small.view { v -> "$v is small" }
result.large.view { v -> "$v is large" }
result.other.view { v -> "$v is other" }