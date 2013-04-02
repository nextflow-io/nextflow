import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel

/**
 * --== ChainWith ==--
 *
 *  A more generic alternative to 'pipe operator
 */

def toUpperCase = {s -> s.toUpperCase()}
final encrypt = new DataflowQueue()
final DataflowReadChannel encrypted = encrypt.chainWith toUpperCase chainWith {it.reverse()} chainWith {'###encrypted###' + it + '###'}

encrypt << "I need to keep this message secret!"
encrypt << "GPars can build linear operator pipelines really easily"

println encrypted.val
println encrypted.val