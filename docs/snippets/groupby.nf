channel.of( tuple(1, 'A'), tuple(1, 'B'), tuple(2, 'C'), tuple(3, 'B'), tuple(1, 'C'), tuple(2, 'A'), tuple(3, 'D') )
    .groupBy()
    .map { key, values -> tuple(key, values.toSorted()) }
    .view()
