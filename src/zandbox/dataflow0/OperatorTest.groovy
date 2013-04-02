import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowVariable



def a = new DataflowVariable()
def b = new DataflowVariable()
def c = new DataflowVariable()
def d = new DataflowVariable()


Dataflow.operator(inputs: [a, b, c], outputs: [d], maxForks: 3) {x, y, z ->
    println "before bind"
    bindOutput 0, x + y + z
    println "after bind"
}


a << 1
b << 2
c << 3

println "before print"
println "=> ${d.val}"

println("Done")