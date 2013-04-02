import groovyx.gpars.dataflow.DataflowQueue

import static groovyx.gpars.dataflow.Dataflow.splitter
import static groovyx.gpars.dataflow.Dataflow.selector

/**
 * Shows how to build operators using the ProcessingNode class
 */

final DataflowQueue aValues = new DataflowQueue()
final DataflowQueue bValues = new DataflowQueue()
final DataflowQueue queue = new DataflowQueue()

splitter(queue, [aValues, bValues])

selector(inputs:[aValues,bValues],outputs:[]) { x, index ->
    println "channel $index => ${x}"
}

queue << 10
queue << 20
queue << 30
queue << 40

sleep 1000