import nextflow.Channel

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

Channel
        .from( 1, 2, 3, 4, 5 )
        .map { it * it  }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }


Channel
        .from( 'a', 'b', 'c' )
        .mapWithIndex { it, index ->  [it, index] }
        .subscribe onNext: { println it }, onComplete: { println 'Done' }

sleep 100