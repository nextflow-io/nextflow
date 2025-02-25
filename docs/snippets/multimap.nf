Channel.of( 1, 2, 3, 4 )
    .multiMap { v ->
        foo: v + 1
        bar: v * v
    }
    .set { result }

result.foo.view { v -> "foo $v" }
result.bar.view { v -> "bar $v" }