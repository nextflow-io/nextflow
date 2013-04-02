import static groovyx.gpars.dataflow.Dataflow.operator

import groovyx.gpars.dataflow.DataflowVariable

/**
 * Shows how to build operators using the ProcessingNode class
 */

final aValues = new DataflowVariable()

//Now the operator is running and processing the data
aValues << 10


operator(inputs: [aValues], outputs:[] )  { x ->
    println x
}.join()
