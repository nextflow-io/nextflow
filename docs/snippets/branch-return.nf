Channel
    .of(1, 2, 3, 40, 50)
    .branch {
        foo: it < 10
            return it+2

        bar: it < 50
            return it-2

        other: true
            return 0
    }