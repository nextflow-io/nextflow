import javax.xml.crypto.Data

import groovyx.gpars.dataflow.DataflowQueue

/**
 * --== Merging channels ==--
 *
 * + Merging allows you to join multiple read channels as inputs for a single dataflow operator.
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

def maleChannel = new DataflowQueue()
def femaleChannel = new DataflowQueue()
def mortgageCandidatesChannel = new DataflowQueue()

maleChannel.merge(femaleChannel) {m, f -> "${m} and ${f}" }.into(mortgageCandidatesChannel)

maleChannel << 'John'
maleChannel << 'Phillip'

femaleChannel << 'Rose'
femaleChannel << 'Kate'



println "Mortgage: " + mortgageCandidatesChannel.getVal()
println "Mortgage: " + mortgageCandidatesChannel.getVal()