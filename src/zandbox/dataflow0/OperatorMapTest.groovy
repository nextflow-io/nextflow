import static groovyx.gpars.dataflow.Dataflow.task

import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.DataflowOperator


/**
 * Shows how to build operators using the ProcessingNode class
 */

aValues = new DataflowQueue()
bValues = new DataflowQueue()
results = new DataflowQueue()


def callback = { println 'ciao' }

def str = '{ x, y -> callback(); bindOutput x + y }'
def closure = new GroovyShell(new Binding(callback:callback)).evaluate(str)

def params = [inputs: [aValues,bValues], outputs: [results], maxForks: 1]
new DataflowOperator(Dataflow.DATA_FLOW_GROUP, params, closure).start()


//Now the operator is running and processing the data
task {

    int count=0
    while(true) {
        aValues << (++count) * 10
        sleep 500
    }

}

task {

    int count=0
    while(true) {
        bValues << (++count)
        sleep 700
    }

}


while(true) println results.val
