import nextflow.Channel

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

numbers = Channel.from( 1, 2, 3 )
results = numbers.mapMany { n -> [ a:n*2, b:n*3 ] }
results.subscribe onNext: { println it }, onComplete: { println 'Done' }

Channel.from ( 1, 2, 3 )
        .mapMany { it -> [ number: it, square: it*it ] }
        .subscribe { println it.key + ': ' + it.value }


sleep 100