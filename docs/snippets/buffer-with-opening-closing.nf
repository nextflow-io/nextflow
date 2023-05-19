// emits bundles starting with `2` and ending with `4`
Channel.of( 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 )
    .buffer( 2, 4 )
    .view()