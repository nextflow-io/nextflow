Channel.of( 1, 2, 3 )
    .multiMap { it -> foo: bar: it }
    .set { result }

result.foo.view { "foo $it" }
result.bar.view { "bar $it" }