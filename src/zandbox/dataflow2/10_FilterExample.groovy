import groovyx.gpars.dataflow.DataflowQueue

/**
 * --== Filtering ==-
 *
 * + The filter() method allows to filter data in the pipeline using boolean predicates.
 *
 */


queue1 = new DataflowQueue()
queue2 = queue1.filter {num -> num % 2 != 0 }

// alternative version:
// queue1.filter(odd) into ( queue2 = new DataflowQueue() )

(1..5).each {queue1 << it}

assert 1 == queue2.val
assert 3 == queue2.val
assert 5 == queue2.val