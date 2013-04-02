import static groovyx.gpars.dataflow.Dataflow.operator

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.PoisonPill

/**
 *  --== Operators ==--
 *
 */

queue = new DataflowQueue<>()
def op = operator( inputs: [queue], outputs: [] ) { println it }

queue << 10 << 20  << PoisonPill.instance


op.join()