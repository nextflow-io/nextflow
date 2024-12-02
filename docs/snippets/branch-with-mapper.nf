Channel.of(1, 2, 3, 40, 50)
    .branch { v ->
        foo: v < 10
            return v + 2

        bar: v < 50
            return v - 2

        other: true
            return 0
    }
    .set { result }

result.foo.view { v -> "$v is foo" }
result.bar.view { v -> "$v is bar" }
result.other.view { v -> "$v is other" }