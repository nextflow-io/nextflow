def criteria = multiMapCriteria { v ->
    small: [v, v < 10]
    large: [v, v > 10]
}

Channel.of(1, 2, 30).multiMap(criteria).set { ch1 }
Channel.of(10, 20, 1).multiMap(criteria).set { ch2 }

ch1.small.view { v, is_small -> "ch1: $v is small: $is_small" }
ch1.large.view { v, is_large -> "ch1: $v is large: $is_large" }
ch2.small.view { v, is_small -> "ch2: $v is small: $is_small" }
ch2.large.view { v, is_large -> "ch2: $v is large: $is_large" }