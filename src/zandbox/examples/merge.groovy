import nextflow.Channel

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */


odds  = Channel.from([1, 3, 5, 7, 9]);
evens = Channel.from([2, 4, 6]);

odds.merge( evens ) { o, e -> [o, e] }
        .subscribe { println it }


sleep 100