Channel.of(
        [1, [1], ['A']],
        [2, [1, 2], ['B', 'C']],
        [3, [1, 2, 3], ['D', 'E']]
    )
    .transpose()
    .view()