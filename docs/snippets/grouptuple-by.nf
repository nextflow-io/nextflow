Channel.of( [1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'] )
    .groupTuple(by: 1)
    .view()