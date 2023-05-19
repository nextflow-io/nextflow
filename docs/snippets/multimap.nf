Channel.of( 1, 2, 3, 4 )
    .multiMap { it ->
        foo: it + 1
        bar: it * it
    }
    .set { result }

result.foo.view { "foo $it" }
result.bar.view { "bar $it" }