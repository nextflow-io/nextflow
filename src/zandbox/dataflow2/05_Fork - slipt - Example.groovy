import static groovyx.gpars.dataflow.Dataflow.task

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel

/**
 *  --== Forking the data flow ==--
 *
 *  + When a need comes to copy the output of a pipeline/channel into more than one following pipeline/channel,
 *    the split() method will help you:
 *
 */

def toUpperCase = {s -> s.toUpperCase()}

final encrypt = new DataflowQueue()
final DataflowWriteChannel messagesToSave = new DataflowQueue()
final DataflowWriteChannel messagesToLog = new DataflowQueue()
encrypt.chainWith toUpperCase split(messagesToSave, messagesToLog)

encrypt << "I need to keep this message secret!"
encrypt << "GPars can build operator pipelines really easy"


println 'log: ' + messagesToLog.getVal() + ' - '+ messagesToLog.getVal()
println 'save: ' + messagesToSave.getVal() + ' - '+ messagesToSave.getVal()


