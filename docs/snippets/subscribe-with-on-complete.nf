Channel.of( 1, 2, 3 )
    .subscribe onNext: { v -> println v }, onComplete: { println 'Done' }