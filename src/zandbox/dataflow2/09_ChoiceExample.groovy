import groovyx.gpars.dataflow.DataflowQueue


/**
 * --== Choices ==--
 *
 *  The binaryChoice() and choice() methods allow you to send a value to one out of two (or many) output channels,
 *  as indicated by the return value from a closure
 */

queue1 = new DataflowQueue<>()
queue2 = new DataflowQueue<>()
queue3 = new DataflowQueue<>()
queue4 = new DataflowQueue<>()

queue1 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9
queue1.choice([queue2, queue3, queue4]) {a -> a % 3}
//queue1.binaryChoice(queue2, queue3) {a -> a > 0}

println "queue2: ${queue2.getVal()} - ${queue2.getVal()} - ${queue2.getVal()} "
println "queue3: ${queue3.getVal()} - ${queue3.getVal()} - ${queue3.getVal()} "
println "queue4: ${queue4.getVal()} - ${queue4.getVal()} - ${queue4.getVal()} "