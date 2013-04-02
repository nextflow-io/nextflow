

import groovyx.gpars.dataflow.DataflowQueue
import static groovyx.gpars.dataflow.Dataflow.task

def words = ['Groovy', 'fantastic', 'concurrency', 'fun', 'enjoy', 'safe', 'GPars', 'data', 'flow']
buffer = new DataflowQueue()

task {
    for (word in words) {
        sleep 1000
        buffer << word.toUpperCase()  //add to the buffer
    }
}


while(true) {
    print "wait.. "
    println "=> ${buffer.val}"
}
