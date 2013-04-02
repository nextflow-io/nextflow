import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel

/**
 * --== Pipe operator ==--
 *
 * + The 'pipe' operator is an handy shortcuts for the common scenario of building (mostly linear) pipelines
 *
 */
def toUpperCase = {s -> s.toUpperCase()}
final encrypt = new DataflowQueue()
final DataflowReadChannel encrypted = encrypt | toUpperCase | {it.reverse()} | {'###encrypted###' + it + '###'}

encrypt << "I need to keep this message secret!"
encrypt << "GPars can build linear operator pipelines really easily"


println encrypted.val
println encrypted.val