// no condition is specified, emits the very first item: 1
Channel.of( 1, 2, 3 )
    .first()
    .view()

// emits the first item matching the regular expression: 'aa'
Channel.of( 'a', 'aa', 'aaa' )
    .first( ~/aa.*/ )
    .view()

// emits the first String value: 'a'
Channel.of( 1, 2, 'a', 'b', 3 )
    .first( String )
    .view()

// emits the first item for which the predicate evaluates to true: 4
Channel.of( 1, 2, 3, 4, 5 )
    .first { v -> v > 3 }
    .view()