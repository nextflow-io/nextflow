import static groovyx.gpars.dataflow.Dataflow.splitter

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel

def toUpperCase = {s -> s.toUpperCase()}
def save = {text ->
    //Just pretending to be saving the text to disk, database or whatever
    println 'Saving ' + text
}

toEncrypt = new DataflowQueue()
encrypted = toEncrypt.chainWith toUpperCase chainWith {it.reverse()} chainWith {'###encrypted###' + it + '###'}

fork1 = new DataflowQueue()
fork2 = new DataflowQueue()
splitter(encrypted, [fork1, fork2])  //Split the data flow

fork1.chainWith save  //Hook in the save operation

//Hook in a sneaky decryption pipeline
decrypted = fork2.chainWith {it[15..-4]}
    .chainWith {it.reverse()}
    .chainWith {it.toLowerCase()}
    .chainWith {'Groovy leaks! Check out a decrypted secret message: ' + it}

toEncrypt << "I need to keep this message secret!"
toEncrypt << "GPars can build operator pipelines really easy"

println "result> " + decrypted.val
println "result> " + decrypted.val