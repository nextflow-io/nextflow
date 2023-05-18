def criteria = branchCriteria {
    small: it < 10
    large: it > 10
}

Channel.of(1, 2, 30).branch(criteria).set { ch1 }
Channel.of(10, 20, 1).branch(criteria).set { ch2 }