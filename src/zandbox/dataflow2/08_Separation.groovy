import groovyx.gpars.dataflow.DataflowQueue

/**
 * --== Separation ==--
 *
 * + Separation is the opposite operation to merge. The supplied closure returns a list of values,
 *   each of which will be output into an output channel with the corresponding position index.
 *
 */

queue1 = new DataflowQueue<>()
queue2 = new DataflowQueue<>()
queue3 = new DataflowQueue<>()
queue4 = new DataflowQueue<>()

queue1 << 1 << 2 << 3

queue1.separate([queue2, queue3, queue4]) {a -> [a-1, a, a+1]}

//println "queue1: ${queue1.getVal()} - ${queue1.getVal()} - ${queue1.getVal()} "
println "queue2: ${queue2.getVal()} - ${queue2.getVal()} - ${queue2.getVal()} "
println "queue3: ${queue3.getVal()} - ${queue3.getVal()} - ${queue3.getVal()} "
println "queue4: ${queue4.getVal()} - ${queue4.getVal()} - ${queue4.getVal()} "