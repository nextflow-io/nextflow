def criteria = branchCriteria { v ->
    small: v < 10
    large: v > 10
}

Channel.of(1, 2, 30).branch(criteria).set { ch1 }
Channel.of(10, 20, 3).branch(criteria).set { ch2 }

ch1.small.view { v -> "$v is small" }
ch1.large.view { v -> "$v is large" }
ch2.small.view { v -> "$v is small" }
ch2.large.view { v -> "$v is large" }