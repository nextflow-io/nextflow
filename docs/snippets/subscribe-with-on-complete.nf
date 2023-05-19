Channel.of( 1, 2, 3 )
    .subscribe onNext: { println it }, onComplete: { println 'Done' }