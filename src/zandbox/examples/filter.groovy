import nextflow.Channel

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

Channel
        .from( 1, 2, 3, 4, 5 )
        .filter { it % 2 == 1 }
        .subscribe { println it }


sleep 100