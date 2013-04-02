import groovyx.gpars.dataflow.DataflowQueue
import static groovyx.gpars.dataflow.Dataflow.task

/**
 * Shows how to build operators using the ProcessingNode class
 */

aValues = new DataflowQueue()
bValues = new DataflowQueue()
results = new DataflowQueue()
join = new DataflowQueue()


groovyx.gpars.dataflow.Dataflow.selector(inputs: [aValues,bValues], outputs:[join] )

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


while(true) println join.val
