import groovyx.gpars.dataflow.DataflowQueue

/**
 * --==  Tapping into the pipeline ==--
 *
 * + Like split() the tap() method allows you to fork the data flow into multiple channels.
 *   Tapping, however, is slightly more convenient in some scenarios, since it treats one of
 *   the two new forks as the successor of the pipeline.
 */

def queue = new DataflowQueue()
def logChannel = new DataflowQueue()
def printChannel = new DataflowQueue()

queue.chainWith {it * 2}.chainWith{it + 10}.tap(logChannel).into(printChannel)

queue << 1 << 2 << 3

logChannel.each { println "log: $it"  }

//printChannel.each { println "print: $it"  }

sleep 1000