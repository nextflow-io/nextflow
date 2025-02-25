Channel.of( ['homer', 5], ['bart', 2], ['lisa', 10], ['marge', 3], ['maggie', 7] )
    .toSortedList { a, b -> b[1] <=> a[1] }
    .view()