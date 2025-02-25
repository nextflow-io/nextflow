Channel.of( 1, 2, 3 )
    .multiMap { v -> foo: bar: v }
    .set { result }

result.foo.view { v -> "foo $v" }
result.bar.view { v -> "bar $v" }