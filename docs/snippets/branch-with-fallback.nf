Channel.of(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it < 50
        other: true
    }
    .set { result }

result.small.view { "$it is small" }
result.large.view { "$it is large" }
result.other.view { "$it is other" }