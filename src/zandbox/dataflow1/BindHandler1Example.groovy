import groovyx.gpars.dataflow.DataflowVariable


/**
 *  --== Bind handler ==--
 *
 *  + Bind handlers can be registered on all dataflow channels (variables, queues or broadcasts) either
 *    using the >> operator and the then() or the whenBound() methods.
 *
 *  + They will be run once a value is bound to the variable.
 */

def a = new DataflowVariable()
a >> {println "The variable has just been bound to $it"}
a.whenBound {println "Just to confirm that the variable has been really set to $it"}

a << 1

sleep 1000
