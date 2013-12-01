import nextflow.Channel

Channel
        .from( 1, 2, 3, 4, 5 )
        .reduce(100) { a, b ->  a+b }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

sleep 100