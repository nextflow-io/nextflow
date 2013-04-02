import static groovyx.gpars.dataflow.Dataflow.task

import groovyx.gpars.dataflow.DataflowQueue


/**
 *  --== Joining pipelines ==--
 *
 *  + Two pipelines (or channels) can be connected using the into() method
 *  + The output of the encryption pipeline is directly connected to the input of the
 *    saving pipeline (a single channel in out case).
 */
def toUpperCase = {s -> s.toUpperCase()}

encrypt = new DataflowQueue()
messagesToSave = new DataflowQueue()
encrypt.chainWith toUpperCase chainWith {it.reverse()} into messagesToSave

task1 = task {
    encrypt << "I need to keep this message secret!"
    encrypt << "GPars can build operator pipelines really easy"
}

task2 = task {
        println "Saving " + messagesToSave.getVal()
}

[task1, task2] *.join()