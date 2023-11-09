def criteria = branchCriteria {
    small: it < 10
    large: it > 10
}

Channel.of(1, 2, 30).branch(criteria).set { ch1 }
Channel.of(10, 20, 3).branch(criteria).set { ch2 }

ch1.small.view { "$it is small" }
ch1.large.view { "$it is large" }
ch2.small.view { "$it is small" }
ch2.large.view { "$it is large" }