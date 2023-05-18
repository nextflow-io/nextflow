Channel
    .of(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it < 50
        other: true
    }